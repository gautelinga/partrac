typedef size_t Uint;

typedef std::vector<std::pair<std::array<Uint, 2>, double>> EdgesType;
typedef std::vector<std::pair<std::array<Uint, 3>, double>> FacesType;
typedef std::list<Uint> FacesListType;
typedef std::list<Uint> EdgesListType;
typedef std::list<Uint> NodesListType;
typedef std::vector<FacesListType> Edge2FacesType;
typedef std::vector<EdgesListType> Node2EdgesType;
typedef std::vector<std::map<Uint, double>> InteriorAnglesType;
typedef Eigen::Vector3d Vector3d;
typedef Eigen::Matrix3d Matrix3d;
