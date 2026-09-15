// SPDX-License-Identifier: MIT
#include "beautize/core.hpp"
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <locale>
#include <set>
#include <sstream>

using namespace beautize;
namespace fs=std::filesystem;
namespace {
const char* help=R"(beautize 0.1.0 - constrained GFN-FF cleanup of Gaussian Cartesian inputs

Usage: beautize input.gjf [options]

The GJF connectivity is mandatory and authoritative; no bond guessing.
GJF -1 freeze flags and Opt=ModRedundant (including ModR) F records are read.
Positive bond orders define edges; their magnitudes do not override GFN-FF.
Atom indices are 1-based. Lengths: Angstrom; angles: degrees.

  -o, --output PATH          Output GJF (default: input.beautized.gjf)
  --freeze 1-12,18,21-25     Freeze Cartesian positions; repeatable
  --mobile 40-58             Freeze the complement of this set
                            Existing freezes still apply inside the mobile set
  --bond 1,2[=1.50]          Fix distance at input value or explicit target
  --angle 1,2,3[=120]        Fix angle; repeatable
  --dihedral 1,2,3,4[=180]   Fix periodic torsion; repeatable
  --constraints PATH        Extra ModRedundant F records; repeatable
  --ignore-gjf-constraints   Intentionally ignore BOTH freeze columns and ModR
  --allow-unconstrained     Explicitly allow optimization without any constraints
  --fmax NUMBER             Projected max atom force, eV/Angstrom (default 0.05)
  --threads INTEGER         GFN-FF threads (default 1; needs OpenMP for >1)
  --max-steps INTEGER       Maximum accepted steps (default 500)
  --max-step NUMBER         Maximum actual atom move, Angstrom (default 0.10)
  --constraint-tol NUMBER   B tolerance in Angstrom, A/D in radians (default 1e-8)
  --memory INTEGER          L-BFGS history size, 1..100 (default 12)
  --report PATH             JSON audit (default: output stem + .json)
  --xyz PATH                Also write final XYZ on convergence
  --trajectory PATH         Accepted-step XYZ trajectory (includes step 0)
  --keep-partial            On failure write output-stem.partial.gjf; exit 2
  --overwrite               Permit replacing output files, never source files
  --dry-run                 Parse, expand, and validate constraints; no GFN-FF
  --verbose                 Include GFN-FF initialization diagnostics
  -q, --quiet               Suppress routine progress output
  -h, --help                Show this help
  --version                 Print version

Default: refuse unconstrained runs and refuse existing outputs. Original GJF
text is preserved except changed XYZ tokens. CLI constraints are recorded in
JSON, NOT inserted into the Gaussian route or ModRedundant block.
Exit codes: 0 success; 1 input/runtime/I/O error; 2 optimizer not converged.
)";
struct Options {
    std::string input,output,report,xyz,trajectory,mobile;
    std::vector<std::string> freezes,constraint_files;
    std::vector<std::pair<Kind,std::string>> coordinates;
    OptimizeOptions optimizer;
    double tolerance{1e-8};
    bool ignore{},allow_unconstrained{},partial{},overwrite{},dry{},verbose{},quiet{};
};
Options parse_cli(int argc,char** argv) {
    Options o;
    for(int i=1;i<argc;++i) {
        std::string arg=argv[i];
        std::optional<std::string> attached;
        if(arg.rfind("--",0)==0 && arg.find('=')!=std::string::npos) {
            auto equal=arg.find('=');attached=arg.substr(equal+1);arg.resize(equal);
        }
        auto value=[&](){
            if(attached) {auto text=*attached;attached.reset();return text;}
            if(++i>=argc)throw Error("Missing value for "+arg);
            return std::string(argv[i]);
        };
        if(arg=="-o"||arg=="--output")o.output=value();
        else if(arg=="--report")o.report=value();
        else if(arg=="--xyz")o.xyz=value();
        else if(arg=="--trajectory")o.trajectory=value();
        else if(arg=="--freeze")o.freezes.push_back(value());
        else if(arg=="--mobile") {if(!o.mobile.empty())throw Error("Use --mobile only once");o.mobile=value();}
        else if(arg=="--bond")o.coordinates.emplace_back(Kind::Bond,value());
        else if(arg=="--angle")o.coordinates.emplace_back(Kind::Angle,value());
        else if(arg=="--dihedral")o.coordinates.emplace_back(Kind::Dihedral,value());
        else if(arg=="--constraints")o.constraint_files.push_back(value());
        else if(arg=="--fmax")o.optimizer.fmax=parse_real(value(),arg);
        else if(arg=="--max-step")o.optimizer.max_step=parse_real(value(),arg);
        else if(arg=="--threads")o.optimizer.threads=parse_int(value(),arg);
        else if(arg=="--max-steps")o.optimizer.max_steps=parse_int(value(),arg);
        else if(arg=="--memory")o.optimizer.memory=parse_int(value(),arg);
        else if(arg=="--constraint-tol")o.tolerance=parse_real(value(),arg);
        else if(arg=="--ignore-gjf-constraints")o.ignore=true;
        else if(arg=="--allow-unconstrained")o.allow_unconstrained=true;
        else if(arg=="--keep-partial")o.partial=true;
        else if(arg=="--overwrite")o.overwrite=true;
        else if(arg=="--dry-run")o.dry=true;
        else if(arg=="--verbose")o.verbose=true;
        else if(arg=="--quiet"||arg=="-q")o.quiet=true;
        else if(!arg.empty() && arg[0]=='-')throw Error("Unknown option: "+arg);
        else {if(!o.input.empty())throw Error("Expected a single input GJF");o.input=arg;}
        if(attached)throw Error("Option does not accept a value: "+arg);
    }
    if(o.input.empty())throw Error("Missing input file; run beautize --help");
    if(o.output.empty()) {fs::path p(o.input);o.output=(p.parent_path()/(p.stem().string()+".beautized.gjf")).string();}
    if(o.report.empty()) {fs::path p(o.output);p.replace_extension(".json");o.report=p.string();}
    if(o.tolerance<1e-12 || o.tolerance>1e-2)throw Error("--constraint-tol must be in [1e-12,1e-2]");
    if(o.optimizer.threads<1 || o.optimizer.fmax<=0 || o.optimizer.max_step<=0 || o.optimizer.max_steps<0 || o.optimizer.memory<1 || o.optimizer.memory>100)
        throw Error("Invalid optimization parameters");
    return o;
}
std::string partial_name(const std::string& output) {
    fs::path p(output);return (p.parent_path()/(p.stem().string()+".partial.gjf")).string();
}
void check_output_paths(const Options& o) {
    std::vector<std::string> sources=o.constraint_files;sources.push_back(o.input);
    std::set<fs::path> paths;
    auto check=[&](const std::string& value) {
        if(value.empty())return;
        auto p=fs::weakly_canonical(value);
        for(const auto& s:sources) {
            if(p==fs::weakly_canonical(s) || (fs::exists(p) && fs::exists(s) && fs::equivalent(p,s)))
                throw Error("Refusing to overwrite a source file: "+value);
        }
        if(!paths.insert(p).second)throw Error("Output paths must be distinct: "+value);
        if(fs::exists(p) && !o.overwrite)throw Error("Output already exists (use --overwrite): "+value);
        if(fs::exists(p) && !fs::is_regular_file(p))throw Error("Output is not a regular file: "+value);
        if(p.has_parent_path() && !fs::is_directory(p.parent_path()))throw Error("Output directory does not exist: "+value);
    };
    check(o.output);check(o.report);check(o.xyz);check(o.trajectory);if(o.partial)check(partial_name(o.output));
}
}
int main(int argc,char** argv) {
    std::locale::global(std::locale::classic());
    for(int i=1;i<argc;++i) {
        std::string arg=argv[i];
        if(arg=="--help"||arg=="-h") {std::cout<<help;return 0;}
        if(arg=="--version") {std::cout<<"beautize 0.1.0\n";return 0;}
    }
    try {
        auto o=parse_cli(argc,argv);
        auto doc=parse_gjf(read_text(o.input));
        ConstraintSet constraints(doc,o.ignore);constraints.tolerance=o.tolerance;
        if(!o.mobile.empty()) {
            auto selection=parse_selection(o.mobile,static_cast<int>(doc.atoms.size()));
            std::vector<bool> mobile(doc.atoms.size(),false);for(int a:selection)mobile[a]=true;
            for(std::size_t i=0;i<mobile.size();++i)if(!mobile[i])constraints.frozen[i]=true;
        }
        for(const auto& spec:o.freezes)for(int a:parse_selection(spec,static_cast<int>(doc.atoms.size())))constraints.frozen[a]=true;
        for(const auto& [kind,spec]:o.coordinates)constraints.add_cli(kind,spec);
        for(const auto& path:o.constraint_files) {
            std::istringstream in(read_text(path));std::string line;int number=0;
            while(std::getline(in,line))constraints.add_modredundant(line,doc,path+":"+std::to_string(++number));
        }
        if(o.ignore)std::cerr<<"Warning: GJF freeze columns and ModRedundant constraints were explicitly ignored.\n";
        if(doc.multiplicity!=1)std::cerr<<"Note: multiplicity is preserved in the GJF, but GFN-FF is not a spin-state calculation.\n";
        if(!o.quiet) {
            std::cout<<"beautize 0.1.0 | "<<doc.atoms.size()<<" atoms | "<<doc.bonds.size()<<" explicit bonds\n"
                     <<constraints.frozen_count()<<" Cartesian-frozen atoms, "<<constraints.internals.size()<<" internal constraints\n";
            for(const auto& c:constraints.internals)std::cout<<"  "<<constraint_name(c)<<" = "<<std::setprecision(12)<<c.target*(c.kind==Kind::Bond?1:180/pi)
                <<(c.kind==Kind::Bond?" Angstrom":" degree")<<"  ["<<c.source<<"]\n";
        }
        if(constraints.empty() && !o.allow_unconstrained && !o.dry)
            throw Error("No constraints found. Specify --freeze/--mobile/internal constraints, or explicitly use --allow-unconstrained");
        Vector x=doc.positions;constraints.restore(x);
        Vector change=x;for(std::size_t i=0;i<x.size();++i)change[i]-=doc.positions[i];
        double projection_displacement=max_atom_norm(change);
        if(o.dry) {
            int rank{};constraints.project(x,Vector(x.size(),0),&rank);
            if(!o.quiet)std::cout<<"Dry run passed. Constraint rank: "<<rank<<"; initial projection max move: "<<projection_displacement<<" Angstrom.\n"
                <<"No GFN-FF evaluation or output files were created.\n";
            if(constraints.empty()&&!o.allow_unconstrained)std::cerr<<"Note: an actual run would require --allow-unconstrained.\n";
            return 0;
        }
        check_output_paths(o);
        Gfnff::set_threads(o.optimizer.threads);
        Gfnff calculator(doc,x,o.verbose?2:0);
        if(!o.quiet)std::cout<<"GFN-FF actual neighbor graph verified against the GJF.\n"
            <<" step          energy/Eh    projected fmax/eV/A    constraint error\n";
        std::string trajectory;
        auto observer=[&](const StepInfo& h,const Vector& pos) {
            if(!o.quiet)std::cout<<std::setw(5)<<h.step<<" "<<std::fixed<<std::setw(20)<<std::setprecision(10)<<h.energy/hartree_to_ev
                <<" "<<std::setw(20)<<std::setprecision(7)<<h.fmax<<" "<<std::scientific<<std::setprecision(3)<<h.constraint_error<<'\n';
            if(!o.trajectory.empty()) {
                std::ostringstream comment;comment<<std::setprecision(17)<<"beautize step="<<h.step<<" energy_hartree="<<h.energy/hartree_to_ev
                    <<" projected_fmax_ev_per_A="<<h.fmax;
                trajectory+=xyz_frame(doc,pos,comment.str());
            }
        };
        auto result=optimize(x,constraints,[&](const Vector& p){return calculator(p);},o.optimizer,observer);
        // Final audit before any structure is published.
        if(constraints.max_error(result.positions)>constraints.tolerance)throw Error("Final constraint audit failed; refusing to publish structure");
        auto output_text=doc.render(result.positions);
        auto roundtrip=parse_gjf(output_text);
        if(constraints.max_error(roundtrip.positions)>constraints.tolerance)throw Error("Output precision audit failed; refusing to publish structure");
        std::string written;
        if(result.converged)written=o.output;
        else if(o.partial)written=partial_name(o.output);
        auto report=make_report(doc,constraints,result,o.optimizer,o.input,written,projection_displacement);
        // Each file is published atomically. The group is not a multi-file transaction.
        if(!written.empty())write_atomic(written,output_text,o.overwrite);
        write_atomic(o.report,report,o.overwrite);
        if(!o.trajectory.empty())write_atomic(o.trajectory,trajectory,o.overwrite);
        if(result.converged&&!o.xyz.empty())write_atomic(o.xyz,xyz_frame(doc,result.positions,"beautize converged; see JSON audit for constraints"),o.overwrite);
        if(!result.converged) {
            std::cerr<<"Optimization did not converge: "<<result.reason<<". Exit code 2.\n";
            if(!written.empty())std::cerr<<"Partial structure: "<<written<<'\n';
            std::cerr<<"Audit: "<<o.report<<'\n';return 2;
        }
        if(!o.quiet)std::cout<<"Converged. Structure: "<<written<<"\nAudit: "<<o.report<<'\n';
        return 0;
    } catch(const std::exception& e) {std::cerr<<"beautize: "<<e.what()<<'\n';return 1;}
}
