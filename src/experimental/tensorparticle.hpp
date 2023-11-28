#ifndef __EXP_TENSORPARTICLES_HPP
#define __EXP_TENSORPARTICLES_HPP

#include "particles.hpp"

// Definitions
class TensorParticle {
public:
    TensorParticle() : TensorParticle({0., 0., 0.}) { };
    TensorParticle(const Vector& x) : m_x(x), m_u({0., 0., 0.}) {
        // m_edges.reserve(100);
    };
    Vector &     x() { return m_x; };
    Vector &     u() { return m_u; };
    Matrix &     F() { return m_F; };
    Real &       rho() { return m_rho; };
    Real &       p() { return m_p; };
    Real &       c() { return m_c; };
    Real &       tau() { return m_tau; };
    //Vector &     n() { return m_n; };  // directional vector
    //Real &       w() { return m_w; };  // accumulated scalar elongation
    //Real &       S() { return m_S; };  // stretching quantity
    int &        cell_id() { return m_cell_id; };  // cell id (if applicable)
    Vector       get_x() const { return m_x; };
    Vector       get_u() const { return m_u; };
    Real         get_rho() const { return m_rho; };
    Real         get_p() const { return m_p; };
    Matrix       get_F() const { return m_F; };

    //Uint           get_id() const { return m_id; };
    //void           set_id(Uint id) { m_id=id; };
    inline static std::vector<std::string> scalar_fields() {
        std::vector<std::string> fields = {"rho", "p", "c", "tau"}; // , "w", "S"};
        return fields;
    }
    inline static std::vector<std::string> vector_fields() {
        std::vector<std::string> fields = {}; //"n"};
        return fields;
    }
    Real get_scalar(const std::string& field) const {
        if (field == "rho"){
            return m_rho;
        }
        else if (field == "p"){
            return m_p;
        }
        else if (field == "c"){
            return m_c;
        }
        else if (field == "tau"){
            return m_tau;
        }
        //else if (field == "w"){
        //    return m_w;
        //}
        //else if (field == "S"){
        //    return m_S;
        //}
        // or throw exception?
        return 0.;
    }
    Vector get_vector(const std::string& field) const {
        //if (field == "n"){
        //    return m_n;
        //}
        return {0., 0., 0.};
    }
protected:
    Vector      m_x = {0., 0., 0.};
    Vector      m_u = {0., 0., 0.};
    Matrix      m_F = {0., 0., 0.};
    Real        m_rho = 0.;
    Real        m_p = 0.;
    Real        m_c = 0.;
    Real        m_tau = 0.;
    //Vector      m_n = {0., 0., 0.};
    //Real        m_w = 0.;
    //Real        m_S = 0.;
    Uint        m_id = 0;
    int         m_cell_id = -1;
    std::set<Uint> m_edge_ids;
    //Particles<Particle>* m_ps;
};

#endif