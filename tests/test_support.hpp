// SPDX-License-Identifier: MIT
#pragma once
#include "beautize/core.hpp"
#include <cmath>
#include <functional>
#include <iostream>
#include <sstream>
#include <string>
using namespace beautize;
inline void require(bool okay,const char* expression,const char* file,int line) {
    if(!okay)throw std::runtime_error(std::string(file)+":"+std::to_string(line)+": "+expression);
}
#define CHECK(...) require(static_cast<bool>((__VA_ARGS__)),#__VA_ARGS__,__FILE__,__LINE__)
inline void near(double actual,double expected,double tolerance) {
    if(!std::isfinite(actual)||std::abs(actual-expected)>tolerance) {
        std::ostringstream out;out.precision(17);out<<actual<<" != "<<expected<<" +/- "<<tolerance;
        throw std::runtime_error(out.str());
    }
}
inline void expect_error(const std::function<void()>& f,const std::string& needle="") {
    try {f();} catch(const Error& e) {
        if(!needle.empty() && std::string(e.what()).find(needle)==std::string::npos)
            throw std::runtime_error("Unexpected error: "+std::string(e.what()));
        return;
    }
    throw std::runtime_error("Expected an Error but none was thrown");
}
struct Suite {
    int passed{},failed{};
    void test(const std::string& name,const std::function<void()>& body) {
        try{body();++passed;std::cout<<"PASS "<<name<<'\n';}
        catch(const std::exception& e){++failed;std::cerr<<"FAIL "<<name<<": "<<e.what()<<'\n';}
    }
    int finish()const {std::cout<<passed<<" passed, "<<failed<<" failed\n";return failed?1:0;}
};
inline std::string water_text(bool modr=false) {
    return std::string("%chk=keep.chk\n#p ")+(modr?"opt=(ModR,maxcycles=80) ":"")+"hf/3-21g geom=connectivity\n\nA title\ncontinued title\n\n0 1\n"
        "O -1 0.000 0.000 0.000\nH 0 1.18 0.04 0.02\nH 0 -0.32 1.10 -0.08\n\n1 2 1.0 3 1.0\n2\n3\n\n"
        +(modr?"B 1 2 F\nA 2 1 3 F\n\n":"");
}
inline Document geometry(const std::vector<int>& numbers,const Vector& x,const std::vector<Bond>& bonds) {
    Document d;d.positions=x;d.bonds=bonds;
    for(int n:numbers){Atom a;a.number=n;a.label=element_symbol(n);d.atoms.push_back(a);}
    return d;
}
