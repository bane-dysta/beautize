// SPDX-License-Identifier: MIT
#include "beautize/core.hpp"
#include <algorithm>
#include <cctype>
#include <charconv>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <locale>
#include <map>
#include <random>
#include <regex>
#include <set>
#include <sstream>

namespace beautize {
std::string trim(std::string s) {
    auto notspace=[](unsigned char c){ return !std::isspace(c); };
    auto first=std::find_if(s.begin(),s.end(),notspace);
    auto last=std::find_if(s.rbegin(),s.rend(),notspace).base();
    return first<last ? std::string(first,last) : std::string{};
}
std::vector<std::string> words(const std::string& s) {
    std::istringstream in(s); in.imbue(std::locale::classic());
    std::vector<std::string> out; std::string token;
    while(in>>token) out.push_back(token);
    return out;
}
int parse_int(const std::string& s, const std::string& context) {
    int value=0;
    const char* first=s.data();
    if(!s.empty() && s.front()=='+') ++first;
    auto r=std::from_chars(first,s.data()+s.size(),value);
    if(first==s.data()+s.size() || r.ec!=std::errc{} || r.ptr!=s.data()+s.size())
        throw Error(context+": expected integer, got '"+s+"'");
    return value;
}
double parse_real(std::string s, const std::string& context) {
    for(auto& c:s) if(c=='d'||c=='D') c='E';
    std::istringstream in(s); in.imbue(std::locale::classic());
    double value{}; in>>value;
    if(!in || in.peek()!=std::char_traits<char>::eof() || !std::isfinite(value))
        throw Error(context+": expected finite number, got '"+s+"'");
    return value;
}
static const std::vector<std::string>& elements() {
    static const auto e=words("X H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn "
        "Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd "
        "Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr");
    return e;
}
std::string element_symbol(int z) {
    if(z<1 || z>103) throw Error("GFN-FF supports atomic numbers 1 through 103");
    return elements().at(static_cast<std::size_t>(z));
}
static int atom_number(const std::string& label, const std::string& where) {
    auto name=label.substr(0,label.find('('));
    if(name.empty()) throw Error(where+": missing element");
    if(label.find('(')!=std::string::npos && label.back()!=')') throw Error(where+": unbalanced atom metadata");
    if(std::isdigit(static_cast<unsigned char>(name[0]))) {
        int n=parse_int(name,where); element_symbol(n); return n;
    }
    if(name.size()>2 || !std::all_of(name.begin(),name.end(),[](unsigned char c){return std::isalpha(c);}))
        throw Error(where+": unsupported atom label '"+label+"' (no dummy, ghost or ONIOM atoms)");
    name[0]=static_cast<char>(std::toupper(static_cast<unsigned char>(name[0])));
    if(name.size()==2) name[1]=static_cast<char>(std::tolower(static_cast<unsigned char>(name[1])));
    auto it=std::find(elements().begin()+1,elements().end(),name);
    if(it==elements().end()) throw Error(where+": unknown element '"+name+"'");
    return static_cast<int>(it-elements().begin());
}
// Token spans, rather than a regenerated atom line, let the writer change only XYZ.
static std::vector<std::pair<std::size_t,std::size_t>> spans(const std::string& line) {
    std::vector<std::pair<std::size_t,std::size_t>> out;
    std::size_t start=std::string::npos;
    int depth=0;
    for(std::size_t i=0;i<=line.size();++i) {
        char c=i<line.size()?line[i]:' ';
        bool sep=i==line.size() || (depth==0 && (std::isspace(static_cast<unsigned char>(c)) || c==','));
        if(sep) {
            if(start!=std::string::npos) {out.emplace_back(start,i-start); start=std::string::npos;}
        } else {
            if(start==std::string::npos) start=i;
            if(c=='(') ++depth;
            if(c==')') --depth;
            if(depth<0) throw Error("Unbalanced parentheses in coordinate line");
        }
    }
    if(depth!=0) throw Error("Unbalanced parentheses in coordinate line");
    return out;
}
Document parse_gjf(const std::string& text) {
    Document d;
    for(std::size_t p=0;p<text.size();) {
        auto e=text.find('\n',p);
        if(e==std::string::npos) {d.lines.push_back({text.substr(p),""}); break;}
        bool cr=e>p && text[e-1]=='\r';
        d.lines.push_back({text.substr(p,e-p-(cr?1:0)),cr?"\r\n":"\n"}); p=e+1;
    }
    if(d.lines.empty()) throw Error("Empty GJF file");
    for(const auto& l:d.lines) {
        std::string s=trim(l.text);
        std::transform(s.begin(),s.end(),s.begin(),[](unsigned char c){return static_cast<char>(std::tolower(c));});
        if(s=="--link1--") throw Error("Multi-Link1 input is not supported; split it into individual jobs");
    }
    std::size_t p=0;
    auto skip=[&](){while(p<d.lines.size() && trim(d.lines[p].text).empty()) ++p;};
    skip();
    while(p<d.lines.size() && !trim(d.lines[p].text).empty() && trim(d.lines[p].text)[0]=='%') ++p;
    skip();
    if(p==d.lines.size() || trim(d.lines[p].text).front()!='#') throw Error("Missing Gaussian route section (#...)");
    while(p<d.lines.size() && !trim(d.lines[p].text).empty()) d.route+=d.lines[p++].text+" ";
    std::string lower=d.route;
    std::transform(lower.begin(),lower.end(),lower.begin(),[](unsigned char c){return static_cast<char>(std::tolower(c));});
    std::set<std::string> route_words;
    std::regex word("[a-z][a-z0-9]*");
    for(std::sregex_iterator i(lower.begin(),lower.end(),word),end;i!=end;++i) route_words.insert(i->str());
    for(const auto& unsupported:{"oniom","qst2","qst3","addgic","gic","readopt","readoptimize","rdoptimize","readfreeze","modconnect"})
        if(route_words.count(unsupported)) throw Error(std::string("Unsupported Gaussian route option: ")+unsupported);
    if(route_words.count("units") && (route_words.count("bohr")||route_words.count("au")||route_words.count("a0")))
        throw Error("Input must use Angstrom Cartesian coordinates; Units=Bohr/AU is not supported");
    bool modr=false;
    for(const auto& w:route_words) if(w.size()>=4 && std::string("modredundant").compare(0,w.size(),w)==0) modr=true;
    skip();
    if(p==d.lines.size()) throw Error("Missing title and molecular specification");
    while(p<d.lines.size() && !trim(d.lines[p].text).empty()) ++p; // title, preserved verbatim
    skip();
    if(p==d.lines.size()) throw Error("Missing charge and multiplicity");
    auto cm=words(d.lines[p++].text);
    if(cm.size()!=2) throw Error("Expected one charge/multiplicity pair (fragment charge lists are not supported)");
    d.charge=parse_int(cm[0],"Charge"); d.multiplicity=parse_int(cm[1],"Multiplicity");
    if(d.multiplicity<1) throw Error("Multiplicity must be positive");
    while(p<d.lines.size() && !trim(d.lines[p].text).empty()) {
        const auto& line=d.lines[p].text;
        auto sp=spans(line);
        std::string where="GJF line "+std::to_string(p+1);
        if(sp.size()!=4 && sp.size()!=5) throw Error(where+": expected element [0|-1] x y z; only Cartesian input is supported");
        auto token=[&](std::size_t i){return line.substr(sp[i].first,sp[i].second);};
        Atom a; a.label=token(0); a.number=atom_number(a.label,where); a.line=p;
        int off=sp.size()==5?2:1;
        if(off==2) {
            int flag=parse_int(token(1),where+" freeze column");
            if(flag!=0 && flag!=-1) throw Error(where+": freeze column must be 0 or -1");
            a.frozen=flag==-1;
        }
        for(int k=0;k<3;++k) {a.xyz_spans[k]=sp[off+k]; d.positions.push_back(parse_real(token(off+k),where));}
        d.atoms.push_back(a); ++p;
    }
    const int n=static_cast<int>(d.atoms.size());
    if(n==0) throw Error("No Cartesian atoms found");
    skip();
    std::map<std::pair<int,int>,double> bonds;
    // GaussView emits one row for EVERY atom, including bare rows for terminal atoms.
    for(int row=1;row<=n;++row) {
        if(p==d.lines.size() || trim(d.lines[p].text).empty())
            throw Error("Missing connectivity row "+std::to_string(row)+"; explicit GaussView connectivity is required (no guessing)");
        auto w=words(d.lines[p].text);
        std::string where="Connectivity line "+std::to_string(p+1);
        if(w.size()%2!=1 || parse_int(w[0],where)!=row)
            throw Error(where+": expected row "+std::to_string(row)+" followed by atom/bond-order pairs");
        for(std::size_t k=1;k<w.size();k+=2) {
            int j=parse_int(w[k],where);
            double order=parse_real(w[k+1],where);
            if(j<1||j>n||j==row||order<0) throw Error(where+": invalid atom index, self-bond, or negative bond order");
            auto key=std::make_pair(std::min(row-1,j-1),std::max(row-1,j-1));
            auto [it,inserted]=bonds.emplace(key,order);
            if(!inserted && std::abs(it->second-order)>1e-12) throw Error(where+": inconsistent duplicate bond orders");
        }
        ++p;
    }
    if(p<d.lines.size() && !trim(d.lines[p].text).empty()) throw Error("Expected blank line after connectivity block");
    std::vector<int> degree(n,0);
    for(const auto& [ij,order]:bonds) if(order>0) {
        d.bonds.push_back({ij.first,ij.second,order}); ++degree[ij.first]; ++degree[ij.second];
    }
    if(*std::max_element(degree.begin(),degree.end())>40) throw Error("GFN-FF neighbor capacity exceeded (maximum 40 neighbors per atom)");
    skip();
    if(modr) while(p<d.lines.size() && !trim(d.lines[p].text).empty()) {
        d.modredundant.emplace_back(p+1,d.lines[p].text); ++p;
    }
    return d;
}
Vector Document::bond_matrix() const {
    auto n=atoms.size(); Vector m(n*n,0);
    for(const auto& b:bonds) m[b.i*n+b.j]=m[b.j*n+b.i]=b.order;
    return m;
}
std::string Document::render(const Vector& x) const {
    if(x.size()!=positions.size()) throw Error("Writer coordinate count mismatch");
    auto copy=lines;
    for(std::size_t i=0;i<atoms.size();++i) for(int k=2;k>=0;--k) {
        double v=x[3*i+k];
        if(!std::isfinite(v)) throw Error("Refusing to write non-finite coordinates");
        if(v==positions[3*i+k]) continue; // frozen tokens retain their exact original spelling
        std::ostringstream s; s.imbue(std::locale::classic()); s<<std::setprecision(17)<<v;
        auto [begin,len]=atoms[i].xyz_spans[k]; copy[atoms[i].line].text.replace(begin,len,s.str());
    }
    std::string out;
    for(const auto& l:copy) out+=l.text+l.ending;
    return out;
}
std::string read_text(const std::string& path) {
    std::ifstream in(path,std::ios::binary);
    if(!in) throw Error("Cannot open input: "+path);
    std::ostringstream out; out<<in.rdbuf();
    if(in.bad()) throw Error("Cannot read input: "+path);
    return out.str();
}
void write_atomic(const std::string& path,const std::string& contents,bool overwrite) {
    namespace fs=std::filesystem;
    fs::path dest(path);
    if(!overwrite && fs::exists(dest)) throw Error("Output exists (use --overwrite): "+path);
    if(dest.has_parent_path() && !fs::is_directory(dest.parent_path())) throw Error("Output directory does not exist: "+path);
    std::random_device random;
    fs::path temp=dest;
    temp+=".tmp."+std::to_string(std::chrono::steady_clock::now().time_since_epoch().count())+"."+std::to_string(random());
    try {
        std::ofstream out(temp,std::ios::binary|std::ios::trunc);
        if(!out) throw Error("Cannot create temporary output beside: "+path);
        out.write(contents.data(),static_cast<std::streamsize>(contents.size())); out.close();
        if(!out) throw Error("Failed writing output: "+path);
        if(overwrite) fs::rename(temp,dest);
        else {
            // Hard-link publication cannot overwrite a concurrently created file.
            fs::create_hard_link(temp,dest); fs::remove(temp);
        }
    } catch(const std::exception& e) {
        std::error_code ignored; fs::remove(temp,ignored);
        throw Error("Cannot publish "+path+": "+e.what());
    }
}
std::vector<int> parse_selection(const std::string& s,int nat) {
    std::set<int> result;
    std::istringstream in(s); std::string part;
    if(s.empty()||s.back()==',') throw Error("Empty atom selection");
    while(std::getline(in,part,',')) {
        auto dash=part.find('-');
        int a=parse_int(trim(part.substr(0,dash)),"Atom selection");
        int b=dash==std::string::npos?a:parse_int(trim(part.substr(dash+1)),"Atom selection range");
        if(a<1||b<a||b>nat) throw Error("Invalid atom selection '"+part+"'; indices are 1-based");
        for(int i=a;i<=b;++i) result.insert(i-1);
    }
    return {result.begin(),result.end()};
}
std::string xyz_frame(const Document& d,const Vector& x,const std::string& comment) {
    std::ostringstream out; out.imbue(std::locale::classic()); out<<d.atoms.size()<<'\n'<<comment<<'\n'<<std::setprecision(17);
    for(std::size_t i=0;i<d.atoms.size();++i) out<<element_symbol(d.atoms[i].number)<<" "<<x[3*i]<<" "<<x[3*i+1]<<" "<<x[3*i+2]<<'\n';
    return out.str();
}
} // namespace beautize
