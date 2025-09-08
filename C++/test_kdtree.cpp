namespace gmmreg {

void dump_mem_usage()
{
    FILE* f=fopen("/proc/self/statm","rt");
    if (!f) return;
    char str[300];
    size_t n=fread(str,1,200,f);
    str[n]=0;
    printf("MEM: %s\n",str);
    fclose(f);
}

#if 0
template<typename T>
using PointCloud = vnl_matrix<T>;

// And this is the "dataset to kd-tree" adaptor class:
template <typename Derived>
struct PointCloudAdaptor
{
    // typedef typename Derived::coord_t coord_t;
    typedef typename Derived::element_type coord_t;
    //typedef double coord_t;

    const Derived &obj; //!< A const ref to the data set origin

    /// The constructor that sets the data set source
    PointCloudAdaptor(const Derived &obj_) : obj(obj_) { }

    /// CRTP helper method
    inline const Derived& derived() const { return obj; }

    // Must return the number of data points
    inline size_t kdtree_get_point_count() const {
        // return derived().pts.size();
        return derived().rows();
    }

    // Returns the dim'th component of the idx'th point in the class:
    // Since this is inlined and the "dim" argument is typically an immediate value, the
    //  "if/else's" are actually solved at compile time.
    inline coord_t kdtree_get_pt(const size_t idx, const size_t dim) const
    {
        return derived()(idx, dim);
    }

    // Optional bounding-box computation: return false to default to a standard bbox computation loop.
    //   Return true if the BBOX was already computed by the class and returned in "bb" so it can be avoided to redo it again.
    //   Look at bb.size() to find out the expected dimensionality (e.g. 2 or 3 for point clouds)
    template <class BBOX>
    bool kdtree_get_bbox(BBOX& /*bb*/) const { return false; }

}; // end of PointCloudAdaptor
#endif

template <typename T>
void generateRandomPointCloud(PointCloud<T> &point, const size_t N, const double max_range = 10)
{
    std::cout << "Generating "<< N << " point cloud...";
    point.set_size(N, 3);
    for (size_t i = 0; i < N;i++)
    {
        point(i, 0) = max_range * (rand() % 1000) / T(1000.0);
        point(i, 1) = max_range * (rand() % 1000) / T(1000.0);
        point(i, 2) = max_range * (rand() % 1000) / T(1000.0);
    }
    std::cout << "done\n";
}

template <typename T>
void kdtree_demo(const size_t N)
{
    PointCloud<T> cloud;

    // Generate points:
    generateRandomPointCloud<T>(cloud, N);

    T query_pt[3] = { 0.5, 0.5, 0.5 };
    // T query_pt[2] = { 0.5, 0.5};

    //const MatrixAdaptor<T> adaptor(cloud);

    std::unique_ptr<NanoflannTree<T>> tree(new NanoflannTree<T>(cloud));
    tree->tree.buildIndex();

    //typedef PointCloudAdaptor<PointCloud<T>> PC2KD;
    // const PC2KD  pc2kd(cloud); // The adaptor

    // template <typename T, size_t DIM>
#if 0
    using MyKdTree = KDTreeSingleIndexAdaptor<
            L2_Simple_Adaptor<T, MatrixAdaptor<T>> ,
            MatrixAdaptor<T>,
            3>;

    const MatrixAdaptor<T> adaptor(cloud); // The adaptor
    // MyKdTree<T, DIM> index{};
    std::unique_ptr<MyKdTree> index(new MyKdTree(
                3, adaptor, KDTreeSingleIndexAdaptorParams(10 /* max leaf */) ));
    index->buildIndex();
#endif
    dump_mem_usage();

    //auto index = CreateKdTree<T, 3>(adaptor);

    #if 0
    // construct a kd-tree index:
    typedef KDTreeSingleIndexAdaptor<
        L2_Simple_Adaptor<T, PC2KD > ,
        PC2KD // ,
        // 3 /* dim */
        > my_kd_tree_t;
    #endif


    //my_kd_tree_t   index(3 /*dim*/, pc2kd, KDTreeSingleIndexAdaptorParams(10 /* max leaf */) );
    //index.buildIndex();


    dump_mem_usage();

    // do a knn search
    const size_t num_results = 1;
    size_t ret_index;
    T out_dist_sqr;
    nanoflann::KNNResultSet<T> resultSet(num_results);
    resultSet.init(&ret_index, &out_dist_sqr );
    tree->tree.findNeighbors(resultSet, &query_pt[0], nanoflann::SearchParams(10));
    //index.knnSearch(query, indices, dists, num_results, mrpt_flann::SearchParams(10));

    std::cout << "knnSearch(nn="<<num_results<<"): \n";
    std::cout << "ret_index=" << ret_index << " out_dist_sqr=" << out_dist_sqr << endl;


    // ----------------------------------------------------------------
    // knnSearch():  Perform a search for the N closest points
    // ----------------------------------------------------------------
    {
        size_t num_results = 5;
        std::vector<size_t>   ret_index(num_results);
        std::vector<T> out_dist_sqr(num_results);

        num_results = tree->tree.knnSearch(&query_pt[0], num_results, &ret_index[0], &out_dist_sqr[0]);

        // In case of less points in the tree than requested:
        ret_index.resize(num_results);
        out_dist_sqr.resize(num_results);

        cout << "knnSearch(): num_results=" << num_results << "\n";
        for (size_t i = 0; i < num_results; i++)
            cout << "idx["<< i << "]=" << ret_index[i] << " dist["<< i << "]=" << out_dist_sqr[i] << endl;
        cout << "\n";
    }

    // ----------------------------------------------------------------
    // radiusSearch(): Perform a search for the points within search_radius
    // ----------------------------------------------------------------
#if 1
    {
        const T search_radius = static_cast<T>(0.1);
        std::vector<std::pair<size_t, T> >   ret_matches;

        nanoflann::SearchParams params;
        //params.sorted = false;

        const size_t nMatches = tree->tree.radiusSearch(&query_pt[0], search_radius, ret_matches, params);

        cout << "radiusSearch(): radius=" << search_radius << " -> " << nMatches << " matches\n";
        for (size_t i = 0; i < nMatches; i++)
            cout << "idx["<< i << "]=" << ret_matches[i].first << " dist["<< i << "]=" << ret_matches[i].second << endl;
        cout << "\n";
    }
#endif
}

}  // namespace gmmreg

/*
int main()
{
    // Randomize Seed
    srand((unsigned int)time(NULL));
    gmmreg::kdtree_demo<float>(1000000);
    gmmreg::kdtree_demo<double>(1000000);
    return 0;
}
*/

