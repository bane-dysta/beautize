// SPDX-License-Identifier: MIT
#include "beautize/core.hpp"
#include <algorithm>
#include <cmath>
#include <deque>
#include <iomanip>
#include <limits>
#include <locale>
#include <sstream>

namespace beautize {
namespace {
Vector difference(const Vector& a,const Vector& b) {
    Vector d=a;for(std::size_t i=0;i<d.size();++i)d[i]-=b[i];return d;
}
void validate_evaluation(const Evaluation& e,std::size_t size) {
    if(e.gradient.size()!=size || !std::isfinite(e.energy))throw Error("Invalid energy/gradient returned by calculator");
    for(double v:e.gradient)if(!std::isfinite(v))throw Error("Non-finite gradient returned by calculator");
}
struct Correction {Vector s,y;double rho;};
Vector lbfgs_direction(const Vector& g,const std::deque<Correction>& memory) {
    Vector q=g;std::vector<double> alpha(memory.size());
    for(std::size_t k=memory.size();k-->0;) {
        alpha[k]=memory[k].rho*dot(memory[k].s,q);
        for(std::size_t i=0;i<q.size();++i)q[i]-=alpha[k]*memory[k].y[i];
    }
    double scale=1.0/70.0;
    if(!memory.empty()) {
        const auto& last=memory.back();
        scale=std::clamp(dot(last.s,last.y)/dot(last.y,last.y),1e-6,1.0);
    }
    for(auto& v:q)v*=scale;
    for(std::size_t k=0;k<memory.size();++k) {
        double beta=memory[k].rho*dot(memory[k].y,q);
        for(std::size_t i=0;i<q.size();++i)q[i]+=memory[k].s[i]*(alpha[k]-beta);
    }
    for(auto& v:q)v=-v;
    return q;
}
std::string quote(const std::string& s) {
    std::ostringstream out;out<<'"';
    for(unsigned char c:s) {
        switch(c) {
        case '"':out<<"\\\"";break;case '\\':out<<"\\\\";break;
        case '\n':out<<"\\n";break;case '\r':out<<"\\r";break;case '\t':out<<"\\t";break;
        default:if(c<0x20)out<<"\\u"<<std::hex<<std::setw(4)<<std::setfill('0')<<static_cast<int>(c)<<std::dec;else out<<c;
        }
    }
    out<<'"';return out.str();
}
double min_unbonded_distance(const Document& d,const Vector& x) {
    auto b=d.bond_matrix();auto n=d.atoms.size();double minimum=std::numeric_limits<double>::infinity();
    for(std::size_t i=0;i<n;++i)for(std::size_t j=0;j<i;++j)if(b[i*n+j]==0)
        minimum=std::min(minimum,std::hypot(x[3*i]-x[3*j],std::hypot(x[3*i+1]-x[3*j+1],x[3*i+2]-x[3*j+2])));
    return minimum;
}
} // namespace
OptimizeResult optimize(Vector x,const ConstraintSet& c,const std::function<Evaluation(const Vector&)>& evaluate,
                        const OptimizeOptions& o,const std::function<void(const StepInfo&,const Vector&)>& observer) {
    if(o.max_steps<0 || o.memory<1 || o.memory>100 || !std::isfinite(o.fmax) || o.fmax<=0 || !std::isfinite(o.max_step) || o.max_step<=0)
        throw Error("Invalid optimizer parameters");
    c.restore(x);
    OptimizeResult result;
    auto evaluate_checked=[&](const Vector& pos) {
        ++result.evaluations;
        Evaluation e=evaluate(pos);validate_evaluation(e,pos.size());return e;
    };
    Evaluation e=evaluate_checked(x);Vector g=c.project(x,e.gradient);
    auto record=[&](int step,double step_size) {
        StepInfo info{step,e.energy,max_atom_norm(g),c.max_error(x),step_size};
        result.history.push_back(info);if(observer)observer(info,x);
    };
    record(0,0);
    std::deque<Correction> memory;
    for(int iteration=0;iteration<=o.max_steps;++iteration) {
        if(max_atom_norm(g)<=o.fmax && c.max_error(x)<=c.tolerance) {result.converged=true;result.reason="converged";break;}
        if(iteration==o.max_steps) {result.reason="max_steps";break;}
        bool accepted=false;Vector next;Evaluation next_e;std::string last_error;
        for(int attempt=0;attempt<2 && !accepted;++attempt) {
            Vector p;
            if(attempt==0) p=c.project(x,lbfgs_direction(g,memory));
            else {memory.clear();p=g;for(auto& v:p)v*=-1.0/70.0;}
            double magnitude=max_atom_norm(p);
            if(magnitude==0)continue;
            if(magnitude>o.max_step)for(auto& v:p)v*=o.max_step/magnitude;
            double slope=dot(g,p);
            if(slope>=0)continue;
            for(int ls=0;ls<28;++ls) {
                double alpha=std::ldexp(1.0,-ls);Vector trial=x;
                for(std::size_t i=0;i<x.size();++i)trial[i]+=alpha*p[i];
                try {
                    c.restore(trial);
                    // Retraction may amplify a tangential step: bound actual displacement too.
                    double actual_step=max_atom_norm(difference(trial,x));
                    if(actual_step>o.max_step*(1+1e-10))continue;
                    if(actual_step<1e-13)break;
                    Evaluation trial_e=evaluate_checked(trial);
                    if(trial_e.energy<=e.energy+1e-4*alpha*slope) {
                        next=std::move(trial);next_e=std::move(trial_e);accepted=true;break;
                    }
                } catch(const Error& ex) {last_error=ex.what();}
            }
        }
        if(!accepted) {
            result.reason="line_search_failed";
            if(!last_error.empty())result.reason+=": "+last_error;
            break;
        }
        Vector next_g=c.project(next,next_e.gradient);
        auto displacement=difference(next,x);
        // Store tangent secants. Old ambient pairs are reprojected by the final
        // direction projection; curvature guards restart on ill-behaved updates.
        Vector s=c.project(next,displacement);
        Vector y=difference(next_g,c.project(next,g));
        double sy=dot(s,y),ss=dot(s,s),yy=dot(y,y);
        if(sy>1e-10*std::sqrt(ss*yy) && sy>1e-16 && yy>0) {
            memory.push_back({std::move(s),std::move(y),1/sy});
            if(static_cast<int>(memory.size())>o.memory)memory.pop_front();
        } else memory.clear();
        double step_size=max_atom_norm(displacement);
        x=std::move(next);e=std::move(next_e);g=std::move(next_g);
        record(iteration+1,step_size);
    }
    result.positions=std::move(x);return result;
}
std::string make_report(const Document& d,const ConstraintSet& c,const OptimizeResult& r,const OptimizeOptions& o,
                        const std::string& input,const std::string& output,double projection_displacement) {
    if(r.history.empty())throw Error("Cannot report an optimization without evaluations");
    std::ostringstream out;out.imbue(std::locale::classic());out<<std::setprecision(17);
    double frozen_drift=0,max_displacement=0;
    Vector displacement=difference(r.positions,d.positions);
    max_displacement=max_atom_norm(displacement);
    for(std::size_t i=0;i<c.frozen.size();++i)if(c.frozen[i])for(int k=0;k<3;++k)frozen_drift=std::max(frozen_drift,std::abs(displacement[3*i+k]));
    int rank{};c.project(r.positions,Vector(r.positions.size(),0),&rank);
    out<<"{\n  \"schema_version\": 1,\n  \"beautize_version\": \"0.1.0\",\n"
       <<"  \"input\": "<<quote(input)<<",\n  \"output\": "<<(output.empty()?"null":quote(output))<<",\n"
       <<"  \"converged\": "<<(r.converged?"true":"false")<<",\n  \"reason\": "<<quote(r.reason)<<",\n"
       <<"  \"topology_source\": \"gjf_connectivity\",\n  \"actual_gfnff_graph_verified\": true,\n"
       <<"  \"bond_order_policy\": \"positive entries define edges; GFN-FF derives its own force-field bond orders\",\n"
       <<"  \"output_gjf_constraints_modified\": false,\n"
       <<"  \"atom_count\": "<<d.atoms.size()<<",\n  \"bond_count\": "<<d.bonds.size()<<",\n"
       <<"  \"charge\": "<<d.charge<<",\n  \"multiplicity_preserved_not_used_by_force_field\": "<<d.multiplicity<<",\n"
       <<"  \"frozen_atom_count\": "<<c.frozen_count()<<",\n  \"internal_constraint_count\": "<<c.internals.size()<<",\n"
       <<"  \"internal_constraint_rank_after_cartesian_elimination\": "<<rank<<",\n"
       <<"  \"free_tangent_dimensions\": "<<(3*static_cast<int>(d.atoms.size())-3*c.frozen_count()-rank)<<",\n"
       <<"  \"gfnff_threads_requested\": "<<o.threads<<",\n"
       <<"  \"fmax_tolerance_ev_per_angstrom\": "<<o.fmax<<",\n  \"max_step_angstrom\": "<<o.max_step<<",\n"
       <<"  \"constraint_tolerance_bond_angstrom\": "<<c.tolerance<<",\n  \"constraint_tolerance_angle_radian\": "<<c.tolerance<<",\n"
       <<"  \"initial_constraint_projection_max_displacement_angstrom\": "<<projection_displacement<<",\n"
       <<"  \"initial_feasible_energy_hartree\": "<<r.history.front().energy/hartree_to_ev<<",\n"
       <<"  \"final_energy_hartree\": "<<r.history.back().energy/hartree_to_ev<<",\n"
       <<"  \"final_projected_fmax_ev_per_angstrom\": "<<r.history.back().fmax<<",\n"
       <<"  \"frozen_max_absolute_component_drift_angstrom\": "<<frozen_drift<<",\n"
       <<"  \"max_atom_displacement_angstrom\": "<<max_displacement<<",\n"
       <<"  \"accepted_steps\": "<<r.history.back().step<<",\n  \"energy_evaluations\": "<<r.evaluations<<",\n";
    auto number_or_null=[&](double value){if(std::isfinite(value))out<<value;else out<<"null";};
    out<<"  \"initial_minimum_not_directly_bonded_distance_angstrom\": ";number_or_null(min_unbonded_distance(d,d.positions));
    out<<",\n  \"final_minimum_not_directly_bonded_distance_angstrom\": ";number_or_null(min_unbonded_distance(d,r.positions));
    out<<",\n  \"frozen_atoms_1based\": [";bool first=true;
    for(std::size_t i=0;i<c.frozen.size();++i)if(c.frozen[i]) {if(!first)out<<", ";first=false;out<<i+1;}
    out<<"],\n  \"input_bonds\": [";
    for(std::size_t i=0;i<d.bonds.size();++i) {const auto& b=d.bonds[i];if(i)out<<",";out<<"\n    ["<<b.i+1<<", "<<b.j+1<<", "<<b.order<<"]";}
    out<<"\n  ],\n  \"constraints\": [";
    for(std::size_t i=0;i<c.internals.size();++i) {
        const auto& item=c.internals[i];double scale=item.kind==Kind::Bond?1:180/pi;
        double value=evaluate_coordinate(item.kind,item.atoms,r.positions).value;
        if(i)out<<",";
        out<<"\n    {\"coordinate\": "<<quote(constraint_name(item))<<", \"target\": "<<item.target*scale
           <<", \"final_value\": "<<value*scale<<", \"absolute_error\": "<<std::abs(coordinate_error(item.kind,value,item.target))*scale
           <<", \"unit\": "<<quote(item.kind==Kind::Bond?"angstrom":"degree")<<", \"source\": "<<quote(item.source)<<"}";
    }
    out<<"\n  ],\n  \"history\": [";
    for(std::size_t i=0;i<r.history.size();++i) {
        const auto& h=r.history[i];if(i)out<<",";
        out<<"\n    {\"step\": "<<h.step<<", \"energy_hartree\": "<<h.energy/hartree_to_ev
           <<", \"projected_fmax_ev_per_angstrom\": "<<h.fmax<<", \"constraint_max_error_angstrom_or_radian\": "<<h.constraint_error
           <<", \"actual_step_angstrom\": "<<h.step_size<<"}";
    }
    out<<"\n  ]\n}\n";return out.str();
}
} // namespace beautize
