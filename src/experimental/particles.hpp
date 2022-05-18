#ifndef __EXP_PARTICLES_HPP
#define __EXP_PARTICLES_HPP

//#include "H5Cpp.h"
#include <algorithm>

#include "typedefs.hpp"
#include "io.hpp"

// Declarations
class Particle;
//template<class ParticleType> 
class Edge;
//template<class ParticleType> class Face;
template<class ParticleType> class Particles;

// Definitions
class Particle {
public:
    Particle() : Particle({0., 0., 0.}) { };
    Particle(const Vector& x) : m_x(x), m_u({0., 0., 0.}) {
        // m_edges.reserve(100);
    };
    Vector &     x() { return m_x; };
    Vector &     u() { return m_u; };
    Real &       rho() { return m_rho; };
    Real &       p() { return m_p; };
    Real &       c() { return m_c; };
    Real &       tau() { return m_tau; };
    Vector       get_x() const { return m_x; };
    Vector       get_u() const { return m_u; };
    Real           get_rho() const { return m_rho; };
    Real           get_p() const { return m_p; };

    //Uint           get_id() const { return m_id; };
    //void           set_id(Uint id) { m_id=id; };
    inline static std::vector<std::string> scalar_fields() {
        std::vector<std::string> fields = {"rho", "p", "c", "tau"};
        return fields;
    }
    inline static std::vector<std::string> vector_fields() {
        std::vector<std::string> fields = {"n"};
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
        // or throw exception?
        return 0.;
    }
    Vector get_vector(const std::string& field) const {
        if (field == "n"){
            return m_n;
        }
        return {0., 0., 0.};
    }
    void connect_edge(int i){
        m_edge_ids.insert(i);
    }
    std::set<Uint>& edge_ids() { return m_edge_ids; }
    /*
    void connect(Particles<Particle>* ps){
        m_ps = ps;
    }
    Particles<Particle>* ps() { return m_ps; };
    */
    void replace_edge_ids(std::map<Uint, Uint>& old2new){
        std::set<Uint> old_edge_ids(m_edge_ids.begin(), m_edge_ids.end());
        m_edge_ids.clear();
        for ( auto const & old_edge_id : old_edge_ids ){
            m_edge_ids.insert(old2new[old_edge_id]);
        }
    }
protected:
    Vector      m_x = {0., 0., 0.};
    Vector      m_u = {0., 0., 0.};
    Real        m_rho = 0.;
    Real        m_p = 0.;
    Real        m_c = 0.;
    Real        m_tau = 0.;
    Vector      m_n = {0., 0., 0.};
    Uint        m_id = 0;
    std::set<Uint> m_edge_ids;
    //Particles<Particle>* m_ps;
};

//template<class ParticleType>
class Edge {
public:
    //Edge (ParticleType* a, ParticleType* b, const Real l0) : m_a(a), m_b(b), m_l0(l0) {
    //    m_faces.reserve(2);
    //};
    //Edge(std::vector<ParticleType>& particles, const Uint a_id, const Uint b_id, const Real l0)
    //: m_particles(particles), m_a_id(a_id), m_b_id(b_id), m_l0(l0) {};
    Edge(const Uint a_id, const Uint b_id, const Real l0)
        : m_a_id(a_id), m_b_id(b_id), m_l0(l0) {};
    ~Edge() {
        //m_ps = NULL;
    };
    Real l0() const { return m_l0; };
    //std::set<ParticleType*> particle_ptrs() const {
    //    return {m_a, m_b};
    //};
    std::set<Uint> particle_ids() const {
        return {m_a_id, m_b_id};
    }
    /*
    Uint set_a_id(const Uint a_id){
        Uint a_id_old = m_a_id;
        m_a_id = a_id;
        return a_id_old;
    }
    Uint set_b_id(const Uint b_id){
        Uint b_id_old = m_b_id;
        m_b_id = b_id;
        return b_id_old;
    }
    void set_l0(const Real l0){
        m_l0 = l0;
    }*/
    //Real length() { return (m_a->x()-m_b->x()).norm(); };
    template<typename T>
    Real length(T& ps) const { return vector(ps).norm(); }
    //void set_l0(Real l0) { m_l0 = l0; };
    template<typename T>
    void resize(T& ps, Real ds_max){
        // _Assumes_ every node belongs to exactly one edge! Make a check for this. (assert_is_resizable)
        // conserve local elongation
        // Real rescale_factor = ds / length(ps); // typically < 1
        //if (rescale_factor < 1.0){
        Real ds = length(ps);
        if (ds > ds_max){
            //std::cout << rescale_factor << std::endl;
            // m_l0 *= rescale_factor;
            int n_doublings_loc = ceil(log2(ds / ds_max));
            m_n_doublings += n_doublings_loc;
            //ds /= exp2(n_doublings_loc);
            //Vector x_a = m_a->x();
            //Vector x_b = m_b->x();
            //Vector x0 = 0.5*(x_a + x_b);
            //Vector dx = 0.5*(x_a - x_b);
            //Vector dx = m_a->x() - m_b->x();
            Vector dx = vector(ps);
            //m_a->x() = x0 + dx * rescale_factor;
            //m_b->x() = x0 - dx * rescale_factor;
            //m_b->x() = m_a->x() - dx * rescale_factor; // Ex.: x_a - (x_a - x_b) * 1.0 = x_b
            //ps.particles()[m_b_id].x() = ps.particles()[m_a_id].x() - dx * rescale_factor;
            ps.particles()[m_b_id].x() = ps.particles()[m_a_id].x() - dx / exp2(n_doublings_loc);
            //m_ps->particles()[m_b_id].x() = m_ps->particles()[m_a_id].x() - dx * rescale_factor;
        }
    }
    void connect_face(const Uint a_id){
        m_face_ids.push_back(a_id);
    };
    template<typename T>
    Vector vector(T& ps) const {
        //return (m_a->x() - m_b->x());
        return (ps.particles()[m_a_id].x() - ps.particles()[m_b_id].x());
        //return (m_ps->particles()[m_a_id].x() - m_ps->particles()[m_b_id].x());
    }
    template<typename T>
    Vector midpoint(T& ps) const {
        //return (m_a->x() + m_b->x())/2;
        return (ps.particles()[m_a_id].x() + ps.particles()[m_b_id].x())/2;
        //return (m_ps->particles()[m_a_id].x() + m_ps->particles()[m_b_id].x())/2;
    }
    template<typename T>
    void move(T& ps, const Vector& dx){
        //m_a->x() += dx;
        //m_b->x() += dx;
        ps.particles()[m_a_id].x() += dx;
        ps.particles()[m_b_id].x() += dx;
        //m_ps->particles()[m_a_id].x() += dx;
        //m_ps.particles()[m_b_id].x() += dx;
    }
    template<typename T>
    void split(T& ps){
        Uint c_id = ps.add_particle(midpoint(ps));
        ps.add_edge(c_id, m_b_id, m_l0/2);
        m_b_id = c_id;
        m_l0 /= 2;
        // TODO: split faces
        //for ( auto & face_id : m_face_ids ){
        //    m_particles.face
        //}
    }
    void replace_node_ids(std::map<Uint, Uint>& old2new){
        m_a_id = old2new[m_a_id];
        m_b_id = old2new[m_b_id];
    }
    Uint n_doublings() const {
        return m_n_doublings;
    }
    template<typename T>
    Real logelong(T& ps){
        return log(length(ps)/l0()) + n_doublings() * log(2);
    }
    template<typename T>
    Real elong(T& ps){
        return length(ps)/l0() * exp2(n_doublings());
    }
    //void connect(Particles<ParticleType>* ps){
    //    m_ps = ps;
    //}
    //Particles<ParticleType>* ps() { return m_ps; };
protected:
    //ParticleType*   m_a = nullptr;
    //ParticleType*   m_b = nullptr;
    Uint m_a_id;
    Uint m_b_id;
    //std::vector<ParticleType>& m_particles;
    //std::vector<ParticleType>& m_particles() { return m_ps.particles(); }
    Real              m_l0 = 0.;
    Uint              m_n_doublings = 0;
    //Uint            m_id = NULL;
    //std::vector<Face<ParticleType>*> m_faces;
    std::vector<Uint> m_face_ids;
    //Particles<ParticleType>& m_ps;
};

template<class ParticleType>
class Face {
public:
    //Face(Edge<ParticleType>* a, Edge<ParticleType>* b, Edge<ParticleType>* c, const Real A0) : m_a(a), m_b(b), m_c(c), m_A0(A0) {};
    Face(Particles<ParticleType>& ps, const Uint a_id, const Uint b_id, const Uint c_id, const Real A0)
     : m_ps(ps), m_a_id(a_id), m_b_id(b_id), m_c_id(c_id), m_A0(A0) {};
    ~Face() {
        //delete m_ps;
    }
    Real A0() const { return m_A0; }
    /*std::set<ParticleType*> particle_ptrs() const { 
        std::set<ParticleType*> pptrs;
        pptrs.merge(m_a->particle_ptrs());
        pptrs.merge(m_b->particle_ptrs());
        pptrs.merge(m_c->particle_ptrs());
        return pptrs;
    };*/
    std::set<Uint> particle_ids() const {
        std::set<Uint> pids;
        pids.insert(m_edges()[m_a_id].particle_ids().begin(),
                    m_edges()[m_a_id].particle_ids().end());
        pids.insert(m_edges()[m_b_id].particle_ids().begin(),
                    m_edges()[m_b_id].particle_ids().end());
        pids.insert(m_edges()[m_c_id].particle_ids().begin(),
                    m_edges()[m_c_id].particle_ids().end());
        return pids;
    }
    template<typename T>
    Real area(T& ps){
        //Vector a = m_a->vector();
        //Vector b = m_b->vector();
        Vector a = m_edges()[m_a_id].vector(ps);
        Vector b = m_edges()[m_b_id].vector(ps);
        return a.cross(b).norm()/2;
    }
    template<typename T>
    Real logelong(T& ps){
        return log(area(ps)/A0()) + n_doublings() * log(2);
    }
    Real n_doublings() const {
        return m_n_doublings;
    }
    template<typename T>
    Real elong(T& ps){
        return area(ps)/A0() * exp2(n_doublings());
    }
    //void connect(Particles<ParticleType>* ps){
    //    m_ps = ps;
    //}
    //Particles<ParticleType>* ps() { return m_ps; };
protected:
    //Edge<ParticleType>* m_a = nullptr;
    //Edge<ParticleType>* m_b = nullptr;
    //Edge<ParticleType>* m_c = nullptr;
    Uint m_a_id;
    Uint m_b_id;
    Uint m_c_id;
    Real                  m_A0 = NULL;
    //Uint                m_id = NULL;
    Uint m_n_doublings         = 0;
    //std::vector<Edge<ParticleType>>& m_edges() { return m_ps.edges(); };
    std::vector<Edge>& m_edges() { return m_ps.edges(); };
    Particles<ParticleType>& m_ps;
};

template<class ParticleType>
class Particles {
public:
    Particles(const Uint Nrw_max);
    ~Particles() {};
    Uint add_particle(const Vector& x){
        ParticleType particle(x);
        //particle.connect(this);
        Uint id = m_particles.size();
        m_particles.push_back(particle);
        return id;
    };
    /*void add_edge(ParticleType* a, ParticleType* b, const Real w) {
        Edge edge(a, b, w);
        m_edges.push_back(edge); // copy??
        a->connect_edge(m_edges.size()-1); // (&(m_edges.back()));
        b->connect_edge(m_edges.size()-1);// (&(m_edges.back()));
    };
    void add_face(Edge<ParticleType>* a, Edge<ParticleType>* b, Edge<ParticleType>* c, const Real w) {
        Face face(a, b, c, w);
        m_faces.push_back(face); 
        a->connect_face(&(m_faces.back()));
        b->connect_face(&(m_faces.back()));
        c->connect_face(&(m_faces.back()));
    };*/
    void add_edge(const Uint a_id, const Uint b_id, const Real w) {
        //Edge edge(m_particles, a_id, b_id, w);
        //Edge edge(*this, a_id, b_id, w);
        Edge edge(a_id, b_id, w);
        //edge.connect(this);
        m_edges.push_back(edge); // copy??
        m_particles[a_id].connect_edge(m_edges.size()-1); // (&(m_edges.back()));
        m_particles[b_id].connect_edge(m_edges.size()-1); // (&(m_edges.back()));
    };
    void add_face(const Uint a_id, const Uint b_id, const Uint c_id, const Real w) {
        Face face(*this, a_id, b_id, c_id, w);
        //face.connect(this);
        m_faces.push_back(face); 
        m_edges[a_id].connect_face(m_faces.size()-1);
        m_edges[b_id].connect_face(m_faces.size()-1);
        m_edges[c_id].connect_face(m_faces.size()-1);
    };
    //ParticleType* particle_ptr(Uint i) { return &(m_particles[i]); };
    //Edge<ParticleType>* edge_ptr(Uint i) { return &(m_edges[i]); };
    //Face<ParticleType>* face_ptr(Uint i) { return &(m_faces[i]); };
    int dim() const { return (m_particles.size() <= 0 ? -1 : (m_edges.size() <= 0 ? 0 : (m_faces.size() <= 0 ? 1 : 2))); };
    std::vector<ParticleType>& particles() { return m_particles; };
    //std::vector<Edge<ParticleType>>& edges() { return m_edges; };
    std::vector<Edge>& edges() { return m_edges; };
    std::vector<Face<ParticleType>>& faces() { return m_faces; };
    void dump_hdf5(H5::H5File& h5f, const std::string& groupname, std::map<std::string, bool>& output_fields);
    void color_particles(Real c0, Real c1){
        Uint i = 0;
        Uint n = m_particles.size()-1;
        for ( auto & particle : m_particles ){
            particle.c() = c0 + (c1-c0)*(Real)(i)/n;
            ++i;
        }
    };
    void resize_edges(Real ds){
        for ( auto & edge: m_edges ){
            edge.resize(*this, ds);
        }
    }
    void refine(Real ds_max){
        for ( auto & edge : m_edges ){
            if ( edge.length(*this) > ds_max ){
                edge.split(*this);
                /*
                //Uint b_id_old = m_b_id;
                Uint c_id = add_particle(edge.midpoint());
                Uint b_id = edge.set_b_id(c_id);
                Real l0 = edge.l0()/2;
                edge.set_l0(l0);
                add_edge(c_id, b_id, l0);
                */
            }
        }
    }

    void write_statistics(std::ofstream& statfile, const Real t);
    void remove_edges(const std::set<Uint>& edge_ids){
        std::vector<Uint> sorted_edge_ids(edge_ids.begin(), edge_ids.end()); 
        std::sort(sorted_edge_ids.begin(), sorted_edge_ids.end()); // ascending
        for ( auto edgeit = sorted_edge_ids.rbegin(); 
              edgeit != sorted_edge_ids.rend(); ++edgeit ){ // descending
            //m_edges.erase(m_edges.begin() + *edgeit);
            // TODO: faces!
        }
        m_edges.erase(m_edges.begin());
        std::map<Uint, Uint> new_edge_id;
        for (Uint i=0; i<sorted_edge_ids.size(); ++i){
            new_edge_id[sorted_edge_ids[i]] = i;
        }
        for ( auto & node : m_particles ){
            node.replace_edge_ids(new_edge_id);
        }
    }
    void remove_particles(const std::set<Uint>& particle_ids){
        std::vector<Uint> sorted_particle_ids(particle_ids.begin(), particle_ids.end());
        std::sort(sorted_particle_ids.begin(), sorted_particle_ids.end()); // ascending
        for ( auto nodeit = sorted_particle_ids.rbegin();
              nodeit != sorted_particle_ids.rend(); ++nodeit ){ // descending
            m_particles.erase(m_particles.begin() + *nodeit);
        }

        std::set<Uint> edge_ids;
        for ( auto const & node_id : sorted_particle_ids ){
            std::set<Uint> edge_ids_loc = m_particles[node_id].edge_ids();
            edge_ids.merge(edge_ids_loc);
        }
        remove_edges(edge_ids);

        std::map<Uint, Uint> new_node_id;
        for ( Uint i=0; i<sorted_particle_ids.size(); ++i){
            new_node_id[sorted_particle_ids[i]] = i;
        }
        for ( auto & edge : m_edges ){
            edge.replace_node_ids(new_node_id);
        }
    }
protected:
    std::vector<ParticleType>           m_particles;
    //std::vector<Edge<ParticleType>>     m_edges;
    std::vector<Edge>     m_edges;
    std::vector<Face<ParticleType>>     m_faces;
    Uint                                m_Nrw_max;
    //friend class Interpol;
    //template<class InterpolatorType>
    //friend class Integrator;
    void _number_particles() { 
        Uint i = 0;
        for (auto & particle : m_particles){
            particle.set_id(i);
            ++i;
        }
    };
    void _get_x_data(std::vector<Real>&);
    void _get_u_data(std::vector<Real>&);
    void _get_scalar_data(const std::string&, std::vector<Real>&);
    void _get_vector_data(const std::string&, std::vector<Real>&);
};

template<class ParticleType>
Particles<ParticleType>::Particles(const Uint Nrw_max) : m_Nrw_max(Nrw_max)
{
    std::cout << "Max size: " << m_particles.max_size() << " Nrw_max:" << Nrw_max << std::endl;
    m_particles.reserve(Nrw_max);
}

template<class ParticleType>
void Particles<ParticleType>::dump_hdf5(H5::H5File& h5f, const std::string& groupname, std::map<std::string, bool>& output_fields){
    //_number_particles();
    /*
    if (m_faces.size() > 0){
        hsize_t faces_dims[2];
        faces_dims[0] = m_faces.size();
        faces_dims[1] = 3;
        DataSpace faces_dspace(2, faces_dims);
        std::vector<Uint> faces_arr(faces_dims[0]*faces_dims[1]);
        for ( auto & face : m_faces ){
            for (auto & node : face.particle_ptrs()){
                std::cout << node->get_id() << " ";
            }
            std::cout << std::endl;
        }
    }
    else */ 
    if (m_edges.size() > 0){
        hsize_t edges_dims[2];
        edges_dims[0] = m_edges.size();
        edges_dims[1] = 2;
        H5::DataSpace edges_dspace(2, edges_dims);
        std::vector<Uint> edges_arr;
        edges_arr.reserve(m_edges.size() * 2);
        //(edges_dims[0]*edges_dims[1]);

        std::vector<Real> dl; 
        dl.reserve(m_edges.size());
        // (edges_dims[0]);
        std::vector<Real> dl0; 
        dl0.reserve(m_edges.size());
        //(edges_dims[0]);
        //std::vector<Real> elong;
        //elong.reserve(m_edges.size());
        std::vector<Real> logelong;
        logelong.reserve(m_edges.size());
        std::vector<Uint> doublings;
        doublings.reserve(m_edges.size());

        for (auto & edge : m_edges){
            /*for (auto & node : edge.particle_ptrs() ){
                edges_arr.push_back(node->get_id());
                //std::cout << node->get_id() << " ";
            }*/
            for ( auto & id : edge.particle_ids() ){
                edges_arr.push_back(id);
            }
            //std::cout << std::endl;
            dl.push_back(edge.length(*this));  //dist(edges[iedge].first[0], edges[iedge].first[1]);
            dl0.push_back(edge.l0());  // [iedge] = edges[iedge].second;
            logelong.push_back(edge.logelong(*this));
            doublings.push_back(edge.n_doublings());
        }
        H5::DataSet edges_dset = h5f.createDataSet(groupname + "/edges",
                                               H5::PredType::NATIVE_ULONG,
                                               edges_dspace);
        edges_dset.write(edges_arr.data(), H5::PredType::NATIVE_ULONG);

        scalar_to_h5(h5f, groupname + "/dl", dl);
        scalar_to_h5(h5f, groupname + "/dl0", dl0);
        scalar_to_h5(h5f, groupname + "/logelong", logelong);

        hsize_t doublings_dims[2];
        doublings_dims[0] = doublings.size();
        doublings_dims[1] = 1;
        H5::DataSpace doublings_dspace(2, doublings_dims);
        H5::DataSet doublings_dset = h5f.createDataSet(groupname + "/doublings",
                                                       H5::PredType::NATIVE_ULONG,
                                                       doublings_dspace);
        doublings_dset.write(doublings.data(), H5::PredType::NATIVE_ULONG);
    }
    if (m_particles.size() > 0){
        std::vector<Real> _vtmp;
        _vtmp.reserve(m_particles.size() * 3);
        _get_x_data(_vtmp);
        vector_to_h5(h5f, groupname + "/points", _vtmp, 3);
        if (output_fields["u"]){
            _get_u_data(_vtmp);
            vector_to_h5(h5f, groupname + "/u", _vtmp, 3);
        }
        for ( auto const& field : ParticleType::scalar_fields() ){
            if (output_fields[field]){
                _get_scalar_data(field, _vtmp);
                scalar_to_h5(h5f, groupname + "/" + field, _vtmp);
            }
        }
        for ( auto const& field : ParticleType::vector_fields() ){
            if (output_fields[field]){
                _get_vector_data(field, _vtmp);
                vector_to_h5(h5f, groupname + "/" + field, _vtmp, 3);
            }
        }
    }
}

template<class ParticleType>
void Particles<ParticleType>::_get_x_data(std::vector<Real>& a){
    a.clear();
    for (auto const& particle : m_particles){
        Vector x = particle.get_x();
        a.push_back(x[0]);
        a.push_back(x[1]);
        a.push_back(x[2]);
    }
}

template<class ParticleType>
void Particles<ParticleType>::_get_u_data(std::vector<Real>& a){
    a.clear();
    for (auto const& particle : m_particles){
        Vector u = particle.get_u();
        a.push_back(u[0]);
        a.push_back(u[1]);
        a.push_back(u[2]);
    }
}

template<class ParticleType>
void Particles<ParticleType>::_get_scalar_data(const std::string& field, std::vector<Real>& a){
    a.clear();
    for (auto const& particle : m_particles){
        Real s = particle.get_scalar(field);
        a.push_back(s);
    }
}

template<class ParticleType>
void Particles<ParticleType>::_get_vector_data(const std::string& field, std::vector<Real>& a){
    a.clear();
    for (auto const& particle : m_particles){
        Vector v = particle.get_vector(field);
        a.push_back(v[0]);
        a.push_back(v[1]);
        a.push_back(v[2]);
    }
}

#endif