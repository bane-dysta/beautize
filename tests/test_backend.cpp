// SPDX-License-Identifier: MIT
#include "test_support.hpp"
#include "gfnff_interface_c.h"
#include <algorithm>
#include <limits>
extern "C" void* bt_gfnff_create(int,const int*,const double*,int,const double*,int,int*);
extern "C" void bt_gfnff_destroy(void*);
extern "C" void bt_gfnff_graph(void*,int,int*);

int main() {
    Suite s;
    s.test("real GFN-FF frozen water cleanup",[]{
        auto d=parse_gjf(water_text());ConstraintSet c(d);Gfnff ff(d,d.positions);
        auto r=optimize(d.positions,c,[&](const Vector& x){return ff(x);},{});
        CHECK(r.converged);CHECK(r.history.back().energy<r.history.front().energy);
        for(int k=0;k<3;++k)near(r.positions[k],d.positions[k],0);
        CHECK(ff.graph()==std::vector<int>({0,1,1,1,0,0,1,0,0}));
    });
    s.test("real GFN-FF coupled hard distance and angle",[]{
        auto d=parse_gjf(water_text(true));ConstraintSet c(d);Gfnff ff(d,d.positions);
        auto r=optimize(d.positions,c,[&](const Vector& x){return ff(x);},{});
        CHECK(r.converged);CHECK(c.max_error(r.positions)<=c.tolerance);CHECK(c.max_error(parse_gjf(d.render(r.positions)).positions)<=c.tolerance);
    });
    s.test("GFN-FF gradient sign and Bohr/eV unit conversion",[]{
        auto d=parse_gjf(water_text());Gfnff ff(d,d.positions);auto reference=ff(d.positions);
        for(std::size_t i=0;i<d.positions.size();++i) {
            auto plus=d.positions,minus=d.positions;plus[i]+=1e-5;minus[i]-=1e-5;
            near(reference.gradient[i],(ff(plus).energy-ff(minus).energy)/2e-5,3e-5);
        }
    });
    s.test("stretched declared bond survives distance-based cutoff",[]{
        auto d=parse_gjf(water_text());d.positions[3]=3.0;Gfnff ff(d,d.positions);auto graph=ff.graph();
        CHECK(graph[1]==1);CHECK(graph[2]==1);CHECK(graph[5]==0);
        auto x=d.positions;x[3]=4.0;ff(x);CHECK(ff.graph()==graph);
    });
    s.test("nearby unbonded atoms never acquire spurious edges",[]{
        auto d=geometry({8,1,1,8,1,1},{0,0,0,.96,0,0,-.24,.93,0,0,0,2.4,0,0,.60,.90,0,2.7},{{0,1,1},{0,2,1},{3,4,1},{3,5,1}});
        Gfnff ff(d,d.positions);auto graph=ff.graph();CHECK(graph[4]==0);CHECK(graph[3*6+4]==1);
        ConstraintSet c(d);c.frozen[0]=c.frozen[1]=c.frozen[2]=true;
        OptimizeOptions options;options.max_steps=600;
        auto r=optimize(d.positions,c,[&](const Vector& x){return ff(x);},options);
        CHECK(r.converged);CHECK(r.history.back().energy<r.history.front().energy);CHECK(ff.graph()==graph);
        for(int k=0;k<9;++k)near(r.positions[k],d.positions[k],0);
    });
    s.test("explicit graph reproduces legacy auto topology on clean water",[]{
        auto d=geometry({8,1,1},{0,0,0,.9572,0,0,-.239987,.926627,0},{{0,1,1},{0,2,1}});
        Gfnff ff(d,d.positions);auto e=ff(d.positions);
        int numbers[3]={8,1,1};Vector xyz=d.positions;for(auto& v:xyz)v/=bohr_to_angstrom;
        auto old=c_gfnff_calculator_init(3,numbers,reinterpret_cast<double(*)[3]>(xyz.data()),0,0,"");
        CHECK(old.ptr!=nullptr);double energy{},sigma[3][3];double gradient[3][3];int status{};
        c_gfnff_calculator_singlepoint(&old,3,numbers,reinterpret_cast<double(*)[3]>(xyz.data()),&energy,gradient,sigma,nullptr,&status);
        CHECK(status==0);near(e.energy,energy*hartree_to_ev,1e-10);
        for(int i=0;i<3;++i)for(int k=0;k<3;++k)near(e.gradient[3*i+k],gradient[i][k]*hartree_to_ev/bohr_to_angstrom,1e-9);
        c_gfnff_calculator_deallocate(&old);
    });
    s.test("aromatic bond-order magnitudes preserved but not force-field overrides",[]{
        Vector x;std::vector<int> numbers;std::vector<Bond> bonds;
        for(int i=0;i<6;++i) {double a=i*pi/3;x.insert(x.end(),{1.4*std::cos(a),1.4*std::sin(a),0});numbers.push_back(6);}
        for(int i=0;i<6;++i) {double a=i*pi/3;x.insert(x.end(),{2.48*std::cos(a),2.48*std::sin(a),0});numbers.push_back(1);bonds.push_back({i,(i+1)%6,1.5});bonds.push_back({i,i+6,1});}
        auto aromatic=geometry(numbers,x,bonds);auto plain=aromatic;for(auto& b:plain.bonds)b.order=1;
        Gfnff f1(aromatic,x),f2(plain,x);near(f1(x).energy,f2(x).energy,1e-10);CHECK(f1.graph()==f2.graph());
    });
    s.test("native Fortran initializer rejects invalid exact graphs",[]{
        int numbers[3]={8,1,1};double xyz[9]={0,0,0,1.8,0,0,-.45,1.75,0};
        for(int scenario=0;scenario<4;++scenario) {
            Vector bonds={0,1,1,1,0,0,1,0,0};
            if(scenario==0)bonds[1]=0;
            if(scenario==1)bonds[0]=1;
            if(scenario==2)bonds[1]=bonds[3]=-1;
            if(scenario==3)bonds[1]=std::numeric_limits<double>::quiet_NaN();
            int status{};auto handle=bt_gfnff_create(3,numbers,xyz,0,bonds.data(),0,&status);
            CHECK(status!=0);CHECK(handle==nullptr);
        }
    });
    s.test("coincident atoms fail without invented random displacement",[]{
        auto d=parse_gjf(water_text());d.positions[3]=d.positions[4]=d.positions[5]=0;
        expect_error([&]{Gfnff ff(d,d.positions);},"nearly coincident");
    });
    return s.finish();
}
