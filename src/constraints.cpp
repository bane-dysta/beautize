// SPDX-License-Identifier: MIT
#include "beautize/core.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <limits>
#include <numeric>
#include <set>
#include <sstream>

extern "C" void bt_symmetric_eigen(int n,double* a,double* w,int* status);
namespace beautize {
namespace {
// Forward-mode differentiation, local to at most four atoms (12 Cartesian DOFs).
// This is an independent implementation, not an ASE port or numerical Jacobian.
struct Dual {
    double v{}; std::array<double,12> d{};
    Dual(double value=0):v(value) {}
    static Dual variable(double value,int id) {Dual a(value); a.d[id]=1; return a;}
};
Dual operator+(const Dual& a,const Dual& b) {Dual c(a.v+b.v); for(int i=0;i<12;++i)c.d[i]=a.d[i]+b.d[i];return c;}
Dual operator-(const Dual& a,const Dual& b) {Dual c(a.v-b.v); for(int i=0;i<12;++i)c.d[i]=a.d[i]-b.d[i];return c;}
Dual operator*(const Dual& a,const Dual& b) {Dual c(a.v*b.v); for(int i=0;i<12;++i)c.d[i]=a.d[i]*b.v+a.v*b.d[i];return c;}
Dual operator/(const Dual& a,const Dual& b) {
    if(std::abs(b.v)<1e-20) throw Error("Singular internal coordinate: zero denominator");
    Dual c(a.v/b.v); for(int i=0;i<12;++i)c.d[i]=(a.d[i]-c.v*b.d[i])/b.v;return c;
}
Dual root(const Dual& a) {
    if(a.v<1e-24) throw Error("Singular internal coordinate: coincident atoms or collinear angle/dihedral");
    Dual c(std::sqrt(a.v));for(int i=0;i<12;++i)c.d[i]=a.d[i]/(2*c.v);return c;
}
Dual atan_two(const Dual& y,const Dual& x) {
    double denom=x.v*x.v+y.v*y.v;
    if(denom<1e-24) throw Error("Singular internal coordinate: undefined angle/dihedral");
    Dual c(std::atan2(y.v,x.v));for(int i=0;i<12;++i)c.d[i]=(x.v*y.d[i]-y.v*x.d[i])/denom;return c;
}
using D3=std::array<Dual,3>;
D3 sub(const D3& a,const D3& b) {return {a[0]-b[0],a[1]-b[1],a[2]-b[2]};}
D3 mul(const D3& a,const Dual& t) {return {a[0]*t,a[1]*t,a[2]*t};}
Dual scalar(const D3& a,const D3& b) {return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];}
D3 cross(const D3& a,const D3& b) {return {a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]};}
Dual norm(const D3& a) {return root(scalar(a,a));}
void validate_dimensions(const ConstraintSet& c,const Vector& x) {
    if(x.size()!=c.reference.size() || x.size()!=3*c.frozen.size()) throw Error("Constraint coordinate dimension mismatch");
    for(double v:x) if(!std::isfinite(v)) throw Error("Non-finite position supplied to constraints");
}
struct Linearization {
    std::vector<Vector> rows;
    Vector residual, eigenvectors, eigenvalues;
    int rank{};
    Linearization(const ConstraintSet& c,const Vector& x) {
        const int m=static_cast<int>(c.internals.size());
        if(m>3000) throw Error("More than 3000 internal constraints: use Cartesian --freeze for a large core instead");
        for(const auto& item:c.internals) {
            auto q=evaluate_coordinate(item.kind,item.atoms,x);
            for(std::size_t i=0;i<c.frozen.size();++i) if(c.frozen[i])
                std::fill(q.gradient.begin()+3*i,q.gradient.begin()+3*i+3,0);
            double scale=std::sqrt(dot(q.gradient,q.gradient));
            if(scale<1e-12) scale=1; // zero row, e.g. all atoms of this coordinate are fixed
            for(auto& value:q.gradient) value/=scale;
            rows.push_back(std::move(q.gradient));
            residual.push_back(coordinate_error(item.kind,q.value,item.target)/scale);
        }
        eigenvectors.assign(static_cast<std::size_t>(m)*m,0); eigenvalues.resize(m);
        for(int j=0;j<m;++j) for(int i=0;i<=j;++i)
            eigenvectors[i+static_cast<std::size_t>(j)*m]=eigenvectors[j+static_cast<std::size_t>(i)*m]=dot(rows[i],rows[j]);
        if(m) {
            int status{}; bt_symmetric_eigen(m,eigenvectors.data(),eigenvalues.data(),&status);
            if(status) throw Error("LAPACK failed to factor the constraint Gram matrix (code "+std::to_string(status)+")");
            double cutoff=std::max(1e-14,eigenvalues.back()*1e-12);
            for(auto& l:eigenvalues) {if(l>cutoff) ++rank;else l=0;}
        }
    }
    Vector solve(const Vector& rhs) const {
        const std::size_t m=rows.size(); Vector result(m,0);
        for(std::size_t k=0;k<m;++k) if(eigenvalues[k]>0) {
            double a=0;
            for(std::size_t i=0;i<m;++i) a+=eigenvectors[i+k*m]*rhs[i];
            a/=eigenvalues[k];
            for(std::size_t i=0;i<m;++i) result[i]+=eigenvectors[i+k*m]*a;
        }
        return result;
    }
    Vector transpose_times(const Vector& w,std::size_t dimension) const {
        Vector result(dimension,0);
        for(std::size_t i=0;i<rows.size();++i) for(std::size_t j=0;j<dimension;++j) result[j]+=rows[i][j]*w[i];
        return result;
    }
};
std::pair<double,double> errors(const ConstraintSet& c,const Vector& x) {
    double maxerr=0,square=0;
    for(const auto& item:c.internals) {
        double e=coordinate_error(item.kind,evaluate_coordinate(item.kind,item.atoms,x).value,item.target);
        maxerr=std::max(maxerr,std::abs(e));square+=e*e;
    }
    return {maxerr,square};
}
} // namespace

double dot(const Vector& a,const Vector& b) {
    if(a.size()!=b.size()) throw Error("Vector dimension mismatch");
    return std::inner_product(a.begin(),a.end(),b.begin(),0.0);
}
double max_atom_norm(const Vector& x) {
    if(x.size()%3) throw Error("Vector length must be a multiple of three");
    double largest=0;
    for(std::size_t i=0;i<x.size();i+=3) largest=std::max(largest,std::hypot(x[i],std::hypot(x[i+1],x[i+2])));
    return largest;
}
ConstraintValue evaluate_coordinate(Kind kind,const std::vector<int>& ids,const Vector& x) {
    if(ids.size()!=static_cast<std::size_t>(kind)) throw Error("Wrong number of atoms in internal coordinate");
    std::array<D3,4> r;
    for(std::size_t i=0;i<ids.size();++i) {
        if(ids[i]<0 || static_cast<std::size_t>(3*ids[i]+2)>=x.size()) throw Error("Internal-coordinate atom index out of bounds");
        for(int k=0;k<3;++k) r[i][k]=Dual::variable(x[3*ids[i]+k],static_cast<int>(3*i)+k);
    }
    Dual q;
    if(kind==Kind::Bond) q=norm(sub(r[0],r[1]));
    else if(kind==Kind::Angle) {
        auto u=sub(r[0],r[1]),v=sub(r[2],r[1]);
        auto sine=norm(cross(u,v));
        if(sine.v/(norm(u).v*norm(v).v)<1e-7) throw Error("Near-linear constrained angle is singular; freeze Cartesian atoms instead");
        q=atan_two(sine,scalar(u,v));
    } else {
        auto b0=sub(r[0],r[1]),b1=sub(r[2],r[1]),b2=sub(r[3],r[2]);
        auto axis=mul(b1,Dual(1)/norm(b1));
        auto v=sub(b0,mul(axis,scalar(b0,axis))),w=sub(b2,mul(axis,scalar(b2,axis)));
        if(norm(v).v/norm(b0).v<1e-7 || norm(w).v/norm(b2).v<1e-7)
            throw Error("Near-collinear constrained dihedral is singular; freeze Cartesian atoms instead");
        q=atan_two(scalar(cross(axis,v),w),scalar(v,w));
    }
    ConstraintValue out;out.value=q.v;out.gradient.assign(x.size(),0);
    if(!std::isfinite(q.v)) throw Error("Non-finite internal coordinate");
    for(std::size_t i=0;i<ids.size();++i) for(int k=0;k<3;++k) {
        double v=q.d[3*i+k]; if(!std::isfinite(v))throw Error("Non-finite internal-coordinate derivative");
        out.gradient[3*ids[i]+k]+=v;
    }
    return out;
}
double coordinate_error(Kind kind,double value,double target) {
    return kind==Kind::Dihedral?std::remainder(value-target,2*pi):value-target;
}
std::string constraint_name(const Constraint& c) {
    std::string name=c.kind==Kind::Bond?"B":c.kind==Kind::Angle?"A":"D";
    for(int i:c.atoms) name+=" "+std::to_string(i+1);
    return name;
}
ConstraintSet::ConstraintSet(const Document& doc,bool ignore):reference(doc.positions),frozen(doc.atoms.size(),false) {
    if(!ignore) {
        for(std::size_t i=0;i<doc.atoms.size();++i) frozen[i]=doc.atoms[i].frozen;
        for(const auto& [line,text]:doc.modredundant) add_modredundant(text,doc,"GJF line "+std::to_string(line));
    }
}
void ConstraintSet::add(Kind kind,std::vector<int> ids,std::optional<double> target,const std::string& source) {
    if(ids.size()!=static_cast<std::size_t>(kind)) throw Error(source+": wrong number of atom indices");
    std::set<int> unique(ids.begin(),ids.end());
    if(unique.size()!=ids.size()) throw Error(source+": repeated atoms in internal coordinate");
    for(int i:ids) if(i<0 || static_cast<std::size_t>(i)>=frozen.size()) throw Error(source+": atom index out of range (indices start at 1)");
    auto reversed=ids;std::reverse(reversed.begin(),reversed.end());if(reversed<ids) ids=reversed;
    double value=target?*target:evaluate_coordinate(kind,ids,reference).value;
    if(!std::isfinite(value)) throw Error(source+": non-finite constraint target");
    if(kind==Kind::Bond && value<1e-6) throw Error(source+": constrained bond length must exceed 0.000001 Angstrom");
    if(kind==Kind::Angle && (value<1e-7 || value>pi-1e-7)) throw Error(source+": angle target must be strictly between 0 and 180 degrees, away from singular limits");
    if(kind==Kind::Dihedral) value=std::remainder(value,2*pi);
    for(auto& old:internals) if(old.kind==kind && old.atoms==ids) {
        if(std::abs(coordinate_error(kind,old.target,value))>1e-10) throw Error(source+": conflicting target for "+constraint_name(old)+" (earlier: "+old.source+")");
        old.source+="; "+source;return;
    }
    internals.push_back({kind,std::move(ids),value,source});
}
void ConstraintSet::add_cli(Kind kind,const std::string& expr) {
    auto eq=expr.find('=');
    std::string indices=expr.substr(0,eq);
    std::replace(indices.begin(),indices.end(),',',' ');
    auto w=words(indices);std::vector<int> ids;
    for(const auto& token:w) ids.push_back(parse_int(token,"CLI constraint")-1);
    std::optional<double> value;
    if(eq!=std::string::npos) {
        value=parse_real(expr.substr(eq+1),"CLI constraint target");
        if(kind!=Kind::Bond) *value*=pi/180;
    }
    add(kind,ids,value,"CLI "+expr);
}
void ConstraintSet::add_modredundant(const std::string& text,const Document& doc,const std::string& source) {
    auto w=words(text.substr(0,text.find_first_of("!#")));
    if(w.empty())return;
    for(auto& word:w) for(auto& c:word) c=static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    if(w.back()!="F") throw Error(source+": only ModRedundant F (freeze) is supported; refusing to ignore '"+text+"'");
    w.pop_back();
    int count=0;bool typed=false;
    if(!w.empty() && w[0].size()==1) {
        if(w[0]=="X") count=1;
        if(w[0]=="B") count=2;
        if(w[0]=="A") count=3;
        if(w[0]=="D") count=4;
        if(count) {typed=true;w.erase(w.begin());}
    }
    if(!typed) count=static_cast<int>(w.size());
    if(count<1||count>4 || (w.size()!=static_cast<std::size_t>(count) && w.size()!=static_cast<std::size_t>(count+1)))
        throw Error(source+": expected [X|B|A|D] indices [value] F");
    std::optional<double> value;
    if(w.size()==static_cast<std::size_t>(count+1)) {
        if(count==1)throw Error(source+": X accepts only atom indices and F");
        value=parse_real(w.back(),source+" target");w.pop_back();
        if(count>2)*value*=pi/180;
    }
    std::vector<int> pattern;bool wildcard=false;
    for(const auto& t:w) {
        if(t=="*") {pattern.push_back(-1);wildcard=true;}
        else {
            int i=parse_int(t,source)-1;
            if(i<0 || static_cast<std::size_t>(i)>=frozen.size()) throw Error(source+": atom index out of bounds");
            pattern.push_back(i);
        }
    }
    if(count==1) {
        if(wildcard) std::fill(frozen.begin(),frozen.end(),true);
        else frozen[pattern[0]]=true;
        return;
    }
    Kind kind=static_cast<Kind>(count);
    if(!wildcard) {add(kind,pattern,value,source);return;}
    // Wildcards range over explicit-graph bonds/angles/proper torsions only.
    // Concrete tuples can also specify nonbonded distances or improper torsions.
    std::vector<std::vector<int>> neighbors(frozen.size());
    for(const auto& b:doc.bonds) {neighbors[b.i].push_back(b.j);neighbors[b.j].push_back(b.i);}
    std::set<std::vector<int>> matches;
    std::vector<int> path;
    std::function<void(int)> visit=[&](int depth){
        if(depth==count) {
            auto canon=path,reverse=path;std::reverse(reverse.begin(),reverse.end());
            if(reverse<canon)canon=reverse;
            matches.insert(canon);return;
        }
        auto consider=[&](int a) {
            if(pattern[depth]!=-1 && pattern[depth]!=a)return;
            if(std::find(path.begin(),path.end(),a)!=path.end())return;
            path.push_back(a);visit(depth+1);path.pop_back();
        };
        if(depth==0)for(std::size_t a=0;a<frozen.size();++a)consider(static_cast<int>(a));
        else for(int a:neighbors[path.back()])consider(a);
    };
    visit(0);
    if(matches.empty())throw Error(source+": wildcard matched no coordinates in the supplied graph");
    for(const auto& ids:matches)add(kind,ids,value,source+" (wildcard)");
}
bool ConstraintSet::empty() const {return internals.empty() && frozen_count()==0;}
int ConstraintSet::frozen_count() const {return static_cast<int>(std::count(frozen.begin(),frozen.end(),true));}
void ConstraintSet::enforce_frozen(Vector& x) const {
    validate_dimensions(*this,x);
    for(std::size_t i=0;i<frozen.size();++i) if(frozen[i])for(int k=0;k<3;++k)x[3*i+k]=reference[3*i+k];
}
double ConstraintSet::max_error(const Vector& x) const {
    validate_dimensions(*this,x);
    double err=errors(*this,x).first;
    for(std::size_t i=0;i<frozen.size();++i) if(frozen[i])for(int k=0;k<3;++k)err=std::max(err,std::abs(x[3*i+k]-reference[3*i+k]));
    return err;
}
Vector ConstraintSet::project(const Vector& x,const Vector& v,int* rank) const {
    validate_dimensions(*this,x);
    if(v.size()!=x.size())throw Error("Projection vector dimension mismatch");
    Vector out=v;
    for(std::size_t i=0;i<frozen.size();++i)if(frozen[i])for(int k=0;k<3;++k)out[3*i+k]=0;
    if(internals.empty()) {if(rank)*rank=0;return out;}
    Linearization j(*this,x); if(rank)*rank=j.rank;
    Vector rhs;for(const auto& row:j.rows)rhs.push_back(dot(row,out));
    auto normal=j.transpose_times(j.solve(rhs),out.size());
    for(std::size_t i=0;i<out.size();++i)out[i]-=normal[i];
    return out;
}
void ConstraintSet::restore(Vector& x,int max_iterations) const {
    enforce_frozen(x);
    if(internals.empty())return;
    if(!std::isfinite(tolerance)||tolerance<=0)throw Error("Constraint tolerance must be positive and finite");
    for(int iteration=0;iteration<max_iterations;++iteration) {
        auto [err,merit]=errors(*this,x);
        if(err<=tolerance)return;
        Linearization j(*this,x);
        Vector rhs=j.residual;for(auto& v:rhs)v=-v;
        auto step=j.transpose_times(j.solve(rhs),x.size());
        double maxstep=max_atom_norm(step);
        if(maxstep<1e-14)throw Error("Incompatible constraints: nonzero residual with no correcting free direction (check frozen atoms and targets)");
        if(maxstep>0.2)for(auto& v:step)v*=0.2/maxstep;
        bool accepted=false;
        for(double alpha=1;alpha>1e-7;alpha*=0.5) {
            Vector trial=x;for(std::size_t i=0;i<x.size();++i)trial[i]+=alpha*step[i];
            enforce_frozen(trial);
            try {
                auto [newerr,newmerit]=errors(*this,trial);
                if(newerr<=tolerance || newmerit<merit*(1-1e-4*alpha)) {x=std::move(trial);accepted=true;break;}
            } catch(const Error&) { /* shrink away from a singular trial geometry */ }
        }
        if(!accepted)throw Error("Constraint restoration stalled: contradictory constraints or singular geometry");
    }
    throw Error("Constraint restoration did not converge within "+std::to_string(max_iterations)+" iterations");
}
} // namespace beautize
