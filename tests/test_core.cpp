// SPDX-License-Identifier: MIT
#include "test_support.hpp"
#include <algorithm>
#include <limits>
#include <random>

int main() {
    Suite s;
    s.test("GJF freeze columns, graph, charge and verbatim roundtrip",[]{
        auto text=water_text();auto d=parse_gjf(text);
        CHECK(d.atoms.size()==3);CHECK(d.bonds.size()==2);CHECK(d.atoms[0].frozen);CHECK(!d.atoms[1].frozen);
        CHECK(d.charge==0&&d.multiplicity==1);CHECK(d.render(d.positions)==text);
    });
    s.test("Opt=(ModR,...) and automatic F constraints",[]{
        auto d=parse_gjf(water_text(true));ConstraintSet c(d);
        CHECK(c.frozen_count()==1);CHECK(c.internals.size()==2);near(c.max_error(d.positions),0,0);
        ConstraintSet ignored(d,true);CHECK(ignored.empty());
    });
    s.test("multiline route, CRLF, metadata and D exponent",[]{
        std::string text="#p opt=modredundant\r\n hf/3-21g geom=connectivity\r\n\r\nTest\r\n\r\n0 1\r\n"
            "O(Fragment=1) -1 0.0D+00 0 0\r\nH 0 1D0 0 0\r\nH 0 0 1D0 0\r\n\r\n1 2 1 3 1\r\n2\r\n3\r\n\r\nB 1 2 F\r\n";
        auto d=parse_gjf(text);CHECK(d.atoms[0].number==8);CHECK(d.render(d.positions)==text);ConstraintSet c(d);CHECK(c.internals.size()==1);
        Vector changed=d.positions;changed[3]+=0.01;auto output=d.render(changed);
        CHECK(output.find("O(Fragment=1) -1 0.0D+00 0 0\r\n")!=std::string::npos);
        CHECK(output.find("1 2 1 3 1\r\n2\r\n3\r\n")!=std::string::npos);
    });
    s.test("basis/ECP tail and numeric bond orders are untouched",[]{
        auto text=water_text(true);text+="O H 0\n6-31G(d)\n****\n\nO 0\nCustom-ECP\n\n";
        auto pos=text.find("1 2 1.0 3 1.0");text.replace(pos,13,"1 2 1.5 3 1.0");
        auto d=parse_gjf(text);near(d.bonds[0].order,1.5,0);Vector x=d.positions;x[3]+=.02;
        auto output=d.render(x);CHECK(output.substr(output.find("O H 0"))==text.substr(text.find("O H 0")));
        CHECK(output.find("1 2 1.5 3 1.0")!=std::string::npos);
    });
    s.test("reject absent or truncated explicit connectivity",[]{
        auto text=water_text();text.erase(text.find("1 2 1.0 3 1.0"));expect_error([&]{parse_gjf(text);},"connectivity");
    });
    s.test("reject malformed, self and inconsistent duplicate edges",[]{
        auto text=water_text();auto p=text.find("1 2 1.0 3 1.0");
        auto copy=text;copy.replace(p,13,"1 1 1.0 3 1.0");expect_error([&]{parse_gjf(copy);},"self-bond");
        copy=text;copy.replace(p,13,"1 2 1.0 3");expect_error([&]{parse_gjf(copy);},"expected row");
        copy=text;auto second=copy.find("\n2\n3");copy.replace(second,3,"\n2 1 2.0\n");expect_error([&]{parse_gjf(copy);},"duplicate");
    });
    s.test("reject unsupported coordinate/route forms and NaN",[]{
        auto text=water_text();auto q=text.find("hf/3-21g");
        auto copy=text;copy.insert(q,"units=bohr ");expect_error([&]{parse_gjf(copy);},"Angstrom");
        copy=text;copy.insert(q,"opt=qst2 ");expect_error([&]{parse_gjf(copy);},"qst2");
        copy=text;copy+="--Link1--\n";expect_error([&]{parse_gjf(copy);},"Multi-Link1");
        copy=text;auto p=copy.find("1.18");copy.replace(p,4,"NaN");expect_error([&]{parse_gjf(copy);},"finite");
    });
    s.test("selection ranges and strict bounds",[]{
        CHECK(parse_selection("1-3,5,2",6)==std::vector<int>({0,1,2,4}));
        for(const auto& value:{"0","1-8","4-2","1,,2","1,"})expect_error([&]{parse_selection(value,6);});
    });
    s.test("independent automatic Jacobians match central differences",[]{
        std::mt19937 generator(1789);std::uniform_real_distribution<double> dist(-2,2);
        for(int trial=0;trial<40;++trial) {
            Vector x(12);for(auto& value:x)value=dist(generator);
            for(auto kind:{Kind::Bond,Kind::Angle,Kind::Dihedral}) {
                std::vector<int> ids;for(int i=0;i<static_cast<int>(kind);++i)ids.push_back(i);
                auto q=evaluate_coordinate(kind,ids,x);
                for(std::size_t j=0;j<x.size();++j) {
                    auto plus=x,minus=x;plus[j]+=1e-6;minus[j]-=1e-6;
                    double finite=coordinate_error(kind,evaluate_coordinate(kind,ids,plus).value,evaluate_coordinate(kind,ids,minus).value)/2e-6;
                    near(q.gradient[j],finite,2e-6);
                }
                auto reverse=ids;std::reverse(reverse.begin(),reverse.end());
                near(coordinate_error(kind,evaluate_coordinate(kind,reverse,x).value,q.value),0,1e-12);
            }
        }
    });
    s.test("dihedral branch-cut shortest periodic residual",[]{
        near(coordinate_error(Kind::Dihedral,179*pi/180,-179*pi/180),-2*pi/180,1e-14);
        near(coordinate_error(Kind::Dihedral,-179*pi/180,179*pi/180),2*pi/180,1e-14);
    });
    s.test("singular angle and dihedral fail explicitly",[]{
        Vector x={0,0,0,1,0,0,2,0,0,3,0,0};
        expect_error([&]{evaluate_coordinate(Kind::Angle,{0,1,2},x);},"Singular");
        expect_error([&]{evaluate_coordinate(Kind::Dihedral,{0,1,2,3},x);},"Singular");
    });
    s.test("typed explicit targets, inferred types and X freezing",[]{
        auto d=parse_gjf(water_text());ConstraintSet c(d,true);
        c.add_modredundant("1 F",d,"test");c.add_modredundant("B 1 2 1.00 F",d,"test");
        c.add_modredundant("A 2 1 3 110 F",d,"test");CHECK(c.frozen_count()==1);CHECK(c.internals.size()==2);
        Vector x=d.positions;c.restore(x);CHECK(c.max_error(x)<=c.tolerance);near(x[0],0,0);near(x[1],0,0);near(x[2],0,0);
        near(evaluate_coordinate(Kind::Bond,{0,1},x).value,1,1e-8);
        near(evaluate_coordinate(Kind::Angle,{1,0,2},x).value,110*pi/180,1e-8);
    });
    s.test("graph wildcard expansion and deduplication",[]{
        auto d=parse_gjf(water_text());ConstraintSet c(d,true);
        c.add_modredundant("B * * F",d,"wildcards");CHECK(c.internals.size()==2);
        c.add_modredundant("A * 1 * F",d,"wildcards");CHECK(c.internals.size()==3);
        c.add_modredundant("2 1 F",d,"duplicate");CHECK(c.internals.size()==3);
        expect_error([&]{c.add_modredundant("D * * * * F",d,"wildcards");},"matched no");
        c.add_modredundant("X * F",d,"all atoms");CHECK(c.frozen_count()==3);
    });
    s.test("torsion wildcard expansion on an explicit path",[]{
        auto d=geometry({6,6,6,6},{0,1,0,0,0,0,1,0,0,1,1,1},{{0,1,1},{1,2,1},{2,3,1}});
        ConstraintSet c(d);c.add_modredundant("D * 2 3 * F",d,"torsions");CHECK(c.internals.size()==1);
    });
    s.test("reject non-freeze records and conflicting explicit targets",[]{
        auto d=parse_gjf(water_text());ConstraintSet c(d);
        expect_error([&]{c.add_modredundant("B 1 2 S 4 0.1",d,"scan");},"only ModRedundant F");
        expect_error([&]{c.add_modredundant("A 1 2 3 A",d,"activate");},"only ModRedundant F");
        c.add_cli(Kind::Bond,"1,2=1.0");expect_error([&]{c.add_cli(Kind::Bond,"2,1=1.1");},"conflicting");
        expect_error([&]{c.add_cli(Kind::Bond,"1,1=1.0");},"repeated");
        expect_error([&]{c.add_cli(Kind::Angle,"1,2,3=180");},"angle target");
    });
    s.test("all-frozen incompatible target is not relaxed",[]{
        auto d=parse_gjf(water_text());ConstraintSet c(d);std::fill(c.frozen.begin(),c.frozen.end(),true);
        c.add_cli(Kind::Bond,"1,2=2.0");auto x=d.positions;
        expect_error([&]{c.restore(x);},"Incompatible");CHECK(x==d.positions);
    });
    s.test("coupled distance/angle rank deficiency handled by pseudoinverse",[]{
        auto d=parse_gjf(water_text());ConstraintSet c(d,true);
        c.add(Kind::Bond,{0,1},{},"d1");c.add(Kind::Bond,{0,2},{},"d2");c.add(Kind::Bond,{1,2},{},"d3");
        c.add(Kind::Angle,{1,0,2},{},"redundant angle");
        Vector v(d.positions.size(),0.3);int rank{};auto p=c.project(d.positions,v,&rank);CHECK(rank==3);
        auto pp=c.project(d.positions,p);for(std::size_t i=0;i<p.size();++i)near(p[i],pp[i],1e-11);
        Vector x=d.positions;x[3]+=.1;x[7]+=.04;c.restore(x);CHECK(c.max_error(x)<=1e-8);
    });
    s.test("projection is tangent and frozen coordinates are eliminated",[]{
        auto d=parse_gjf(water_text(true));ConstraintSet c(d);Vector v={3,2,1,1,2,3,5,7,11};
        auto p=c.project(d.positions,v);for(int i=0;i<3;++i)near(p[i],0,0);
        for(const auto& item:c.internals)near(dot(evaluate_coordinate(item.kind,item.atoms,d.positions).gradient,p),0,1e-10);
    });
    s.test("explicit dihedral target crosses periodic branch correctly",[]{
        double phi=179*pi/180;
        auto d=geometry({6,6,6,6},{0,1,0,0,0,0,1,0,0,1,std::cos(phi),std::sin(phi)},{{0,1,1},{1,2,1},{2,3,1}});
        ConstraintSet c(d);c.frozen[0]=c.frozen[1]=c.frozen[2]=true;
        c.add_cli(Kind::Dihedral,"1,2,3,4=-179");Vector x=d.positions;c.restore(x);
        CHECK(c.max_error(x)<=c.tolerance);
        Vector diff=x;for(std::size_t i=0;i<x.size();++i)diff[i]-=d.positions[i];CHECK(max_atom_norm(diff)<0.05);
    });
    s.test("constrained L-BFGS converges and obeys every accepted-step bound",[]{
        auto d=parse_gjf(water_text());ConstraintSet c(d);c.add_cli(Kind::Bond,"1,2=1.0");
        Vector target={0,0,0,0,1.0,0,-.2,.9,0};OptimizeOptions options;options.fmax=1e-6;options.max_steps=400;
        auto potential=[&](const Vector& x){Evaluation e;e.gradient.resize(x.size());for(std::size_t i=0;i<x.size();++i){double q=x[i]-target[i];e.energy+=5*q*q;e.gradient[i]=10*q;}return e;};
        auto r=optimize(d.positions,c,potential,options,[&](const StepInfo& h,const Vector& x){CHECK(c.max_error(x)<=c.tolerance);CHECK(h.step_size<=options.max_step*(1+1e-9));});
        CHECK(r.converged);CHECK(r.history.back().fmax<=options.fmax);
        for(int i=0;i<3;++i)near(r.positions[i],d.positions[i],0);
        for(std::size_t i=1;i<r.history.size();++i)CHECK(r.history[i].energy<=r.history[i-1].energy);
    });
    s.test("nonconvergence is explicit, not silent success",[]{
        auto d=parse_gjf(water_text());ConstraintSet c(d);OptimizeOptions options;options.max_steps=0;
        auto potential=[](const Vector& x){return Evaluation{dot(x,x),x};};
        auto r=optimize(d.positions,c,potential,options);CHECK(!r.converged);CHECK(r.reason=="max_steps");CHECK(r.history.size()==1);
    });
    s.test("fully frozen feasible geometry is a valid constrained solution",[]{
        auto d=parse_gjf(water_text(true));ConstraintSet c(d);std::fill(c.frozen.begin(),c.frozen.end(),true);
        auto r=optimize(d.positions,c,[](const Vector& x){return Evaluation{dot(x,x),Vector(x.size(),123.0)};},{});
        CHECK(r.converged);CHECK(r.positions==d.positions);CHECK(r.history.back().step==0);near(r.history.back().fmax,0,0);
    });
    return s.finish();
}
