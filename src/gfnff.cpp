// SPDX-License-Identifier: MIT
#include "beautize/core.hpp"
#include <cmath>
#include <string>
extern "C" {
void* bt_gfnff_create(int n,const int* at,const double* xyz,int charge,const double* bonds,int level,int* status);
void bt_gfnff_destroy(void* handle);
void bt_gfnff_threads(int n,int* status);
void bt_gfnff_evaluate(void* handle,int n,const int* at,const double* xyz,double* e,double* g,int* status);
void bt_gfnff_graph(void* handle,int n,int* adjacency);
}
namespace beautize {
static Vector to_bohr(const Vector& x) {
    Vector b=x;
    for(auto& v:b) {if(!std::isfinite(v))throw Error("Non-finite GFN-FF input coordinate");v/=bohr_to_angstrom;}
    // GFN-FF cannot define a force direction at coincident atoms. Do not introduce
    // arbitrary random motion into a reference structure without user consent.
    for(std::size_t i=0;i<b.size()/3;++i)for(std::size_t j=0;j<i;++j) {
        double distance=std::hypot(b[3*i]-b[3*j],std::hypot(b[3*i+1]-b[3*j+1],b[3*i+2]-b[3*j+2]));
        if(distance<1e-3)throw Error("Atoms "+std::to_string(j+1)+" and "+std::to_string(i+1)+" are nearly coincident (<0.000529 Angstrom); separate them before optimization");
    }
    return b;
}
Gfnff::Gfnff(const Document& doc,const Vector& x,int printlevel) {
    if(x.size()!=doc.positions.size())throw Error("GFN-FF coordinate dimension mismatch");
    for(const auto& a:doc.atoms)numbers_.push_back(a.number);
    auto xyz=to_bohr(x),bonds=doc.bond_matrix();int status{};
    handle_=bt_gfnff_create(static_cast<int>(numbers_.size()),numbers_.data(),xyz.data(),doc.charge,bonds.data(),printlevel,&status);
    if(status||!handle_)throw Error("GFN-FF initialization failed (code "+std::to_string(status)+"); input graph was NOT replaced with a guessed graph");
    auto actual=graph();
    for(std::size_t i=0;i<bonds.size();++i) if(actual[i]!=(bonds[i]>0?1:0)) {
        bt_gfnff_destroy(handle_);handle_=nullptr;
        throw Error("Internal error: GFN-FF final neighbor list differs from the supplied GJF graph");
    }
}
void Gfnff::set_threads(int count) {
    int status{};bt_gfnff_threads(count,&status);
    if(status)throw Error("Requested thread count is invalid, or this build has no OpenMP (use --threads 1)");
}
Gfnff::~Gfnff(){if(handle_)bt_gfnff_destroy(handle_);}
Evaluation Gfnff::operator()(const Vector& x) {
    if(x.size()!=3*numbers_.size())throw Error("GFN-FF coordinate dimension mismatch");
    auto xyz=to_bohr(x); Evaluation e;e.gradient.resize(x.size());int status{};
    bt_gfnff_evaluate(handle_,static_cast<int>(numbers_.size()),numbers_.data(),xyz.data(),&e.energy,e.gradient.data(),&status);
    if(status)throw Error("GFN-FF energy/gradient failed (code "+std::to_string(status)+")");
    e.energy*=hartree_to_ev;
    for(auto& v:e.gradient)v*=hartree_to_ev/bohr_to_angstrom;
    return e;
}
std::vector<int> Gfnff::graph() const {
    std::vector<int> adjacency(numbers_.size()*numbers_.size());
    bt_gfnff_graph(handle_,static_cast<int>(numbers_.size()),adjacency.data());return adjacency;
}
} // namespace beautize
