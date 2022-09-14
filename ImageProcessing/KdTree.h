#pragma once
// kdtree - An implemention of a kd-tree spatial data structure that
//          is templated on point dimension as well as data type.
//          This class supports building from a vector of points,
//          nearest neighbor queries, serialization to a text file, as
//          well as deserialization.
//
#include <assert.h>
#include <algorithm>
#include <iostream>
#include <memory>
#include <numeric>
#include <sstream>
#include <vector>


// A helper class to track a templated point along with an index that
// refers back to the primal dataset
template<class T>
class IndexedPoint
{
public:
    IndexedPoint() : index(-1)
    {}

    IndexedPoint(int idx, const std::vector<T>& pt)
        : index(idx)
        , point(pt)
    {}

    int index;
    std::vector<T> point;
};


// Forward declarations
//
template<class T> class KDTree;


/////////////////////
// Class declarations
template<class T> class KDNode
{
public:
    KDNode<T>()
        : _leftNode{ nullptr }
        , _rightNode{ nullptr }
    {};

    void build(std::vector< IndexedPoint<T> >& data,
        typename std::vector< IndexedPoint<T> >::iterator start,
        typename std::vector< IndexedPoint<T> >::iterator end,
        int axis = -1);

    void nearestNeighbor(const std::vector<T>& queryPoint,
        IndexedPoint<T>& best, double& bestSqrDistance);

    // We use streams for serialization/deserialization
    template<class U> friend std::ostream& operator<< (std::ostream& out,
        const KDNode<U>& node);
    template<class U> friend std::istream& operator>> (std::istream& is,
        KDNode<U>& node);

    // We want this function that works on trees to be a friend of the
    // node class, so it can access the private constructor below
    template<class U> friend std::istream& operator>> (std::istream& is,
        KDTree<U>& tree);

private:
    // Used during deserialization and made private because we don't
    // want consumers of this class to be able to partially construct
    // nodes in this manner
    KDNode<T>(int axis)
        : _axis(axis)
        , _leftNode{ nullptr }
        , _rightNode{ nullptr }
    {};

    int                         _axis;  // separation axis
    IndexedPoint<T>             _location;

    // pointers, possibly null, to the children of this node
    std::unique_ptr<KDNode<T>>  _leftNode;
    std::unique_ptr<KDNode<T>>  _rightNode;
};


template<class T> class KDTree
{
public:
    KDTree<T>() : _rootNode{ nullptr } {};

    // 'build' takes a vector of points and builds a balanced KD tree.
    bool build(const std::vector< std::vector<T> >& data);

    // Returns the value and index (into the primal dataset) of the
    // closest point (Euclidian distance) to the query point
    IndexedPoint<T> nearestNeighbor(const std::vector<T>& queryPoint);

    // Used in serializaton and deserialization
    template<class U> friend std::ostream& operator<< (std::ostream& out,
        const KDTree<U>& tree);
    template<class U> friend std::istream& operator>> (std::istream& is,
        KDTree<U>& node);

private:
    std::unique_ptr<KDNode<T>> _rootNode;
};



//////////////////
// Implementations

// Helper functions

// Compares two vectors based on their in_axis'th element; used for sorting
template <class T>
class IndexedPointCompare
{
public:
    IndexedPointCompare(int in_axis)
        : _axis{ in_axis }
    {}

    bool operator()(const IndexedPoint<T>& i, const IndexedPoint<T>& j) const
    {
        return i.point[_axis] < j.point[_axis];
    }

private:
    int _axis;
};

// Returns the squared Euclidian distance between two equal-sized
// vectors.
template<class T>
double squaredDistance(const std::vector<T>& a, const std::vector<T>& b)
{
    assert(a.size() == b.size());

    double dist = 0.0;
    for (int ii = 0; ii < a.size(); ++ii)
    {
        dist += (a[ii] - b[ii]) * (a[ii] - b[ii]);
    }

    return dist;
}

// KDTree member function implementations
//
template <class T>
bool KDTree<T>::build(const std::vector< std::vector<T> >& data)
{
    using namespace std;

    bool success = false;
    try
    {
        // Copy the input data over to a non-const vector that also has
        // indices
        vector< IndexedPoint<T> > dataCopy;
        for (int ii = 0; ii < data.size(); ++ii)
        {
            dataCopy.push_back(IndexedPoint<T>(ii, data[ii]));
        }

        _rootNode = std::unique_ptr<KDNode<T>>(new KDNode<T>);
        _rootNode->build(dataCopy, dataCopy.begin(), dataCopy.end());
        success = true;
    }
    catch (exception& e)
    {
        cerr << "Exception during construction: " << e.what() << endl;
    }

    return success;
}

template <class T>
IndexedPoint<T> KDTree<T>::nearestNeighbor(const std::vector<T>& queryPoint)
{
    using namespace std;

    if (_rootNode == nullptr)
    {
        cerr << "No tree has been constructed" << endl;
        return IndexedPoint<T>();
    }

    // Initialize the nearest point to the origin, and assign it a
    // maxiamally bad distance
    IndexedPoint<T> best;
    best.point.resize(queryPoint.size(), 0);
    best.index = -1;
    double dist = std::numeric_limits<T>::max();

    try
    {
        _rootNode->nearestNeighbor(queryPoint, best, dist);
    }
    catch (exception& e)
    {
        cerr << "Exception during nearest neighbor query: " << e.what() << endl;
    }

    return best;
}

// Serialization
template <class T>
std::ostream& operator<< (std::ostream& out, const KDTree<T>& tree)
{
    using namespace std;

    // Ensure that we seialize our results at full numeric precision
    out.precision(numeric_limits<T>::max_digits10);

    if (tree._rootNode != nullptr)
    {
        out << *(tree._rootNode);
    }

    return out;
}

// We use 'axis == -1' as an indictor for a node that doesn't exist.
#define NULLPTR_MARKER -1


static int intFromStream(std::istream& in)
{
    using namespace std;

    // Read an int
    string line;
    getline(in, line);
    istringstream axis_stream(line);

    int val;
    if (!(axis_stream >> val))
    {
        // End of file
        return -1;
    }

    return val;
}

template <class T>
std::ostream& operator<< (std::ostream& out, const KDNode<T>& node)
{
    using namespace std;

    out << node._axis << endl;

    out << node._location.index << endl;
    for (const auto& x : node._location.point)
    {
        out << x << " ";
    }
    out << endl;

    if (node._leftNode != nullptr)
    {
        out << *(node._leftNode);
    }
    else
    {
        out << NULLPTR_MARKER << endl;
    }

    if (node._rightNode != nullptr)
    {
        out << *(node._rightNode);
    }
    else
    {
        out << NULLPTR_MARKER << endl;
    }

    return out;
}


// Deserialization
template<class T>
std::istream& operator>> (std::istream& in, KDTree<T>& tree)
{
    using namespace std;

    // Read the axis so we can determine if our tree should even
    // be constructed
    int axis = intFromStream(in);

    if (axis == -1)
    {
        // A root node doesn't exist or we've reached the end of the file
        return in;
    }

    tree._rootNode = std::unique_ptr<KDNode<T>>(new KDNode<T>(axis));
    in >> *(tree._rootNode);
    return in;
}

template<class T>
std::istream& operator>> (std::istream& in, KDNode<T>& node)
{
    using namespace std;

    node._location.index = intFromStream(in);

    string line;
    getline(in, line);
    istringstream point_stream(line);

    T elem;
    while (point_stream >> elem)
    {
        node._location.point.push_back(elem);
    }

    // read axis of left child
    int axis = intFromStream(in);
    if (axis == -1)
    {
        // The node doesn't exist or we've reached the end of the file
        node._leftNode = std::unique_ptr<KDNode<T>>(nullptr);
    }
    else
    {
        node._leftNode = std::unique_ptr<KDNode<T>>(new KDNode<T>(axis));
        in >> *(node._leftNode);
    }

    // read axis of right child
    axis = intFromStream(in);
    if (axis == -1)
    {
        // The node doesn't exist or we've reached the end of the file
        node._rightNode = std::unique_ptr<KDNode<T>>(nullptr);
    }
    else
    {
        node._rightNode = std::unique_ptr<KDNode<T>>(new KDNode<T>(axis));
        in >> *(node._rightNode);
    }

    return in;
}


// KDNode member function defintion
//
template <class T>
void KDNode<T>::build(std::vector< IndexedPoint<T> >& data,
    typename std::vector< IndexedPoint<T> >::iterator start,
    typename std::vector< IndexedPoint<T> >::iterator end,
    int axis)
{
    using namespace std;

    _axis = (axis + 1) % start->point.size(); // separating axis of this node

    const int dataSize = end - start;
    if (dataSize < 2)
    {
        _location = *start; // copy, so we can discard the data vector
                // after construction
        _leftNode = std::unique_ptr<KDNode<T>>(nullptr);
        _rightNode = std::unique_ptr<KDNode<T>>(nullptr);
        return;
    }

    // find median element
    typename vector<IndexedPoint<T>>::iterator medianItr = start + dataSize / 2;
    nth_element(start, medianItr, end, IndexedPointCompare<T>(_axis));

    _location = *medianItr; // copy, so we can discard the data vector
                // after construction

    typename vector<IndexedPoint<T>>::iterator leftStart = start;
    typename vector<IndexedPoint<T>>::iterator leftEnd = medianItr;
    typename vector<IndexedPoint<T>>::iterator rightStart = medianItr;
    typename vector<IndexedPoint<T>>::iterator rightEnd = end;

    if (leftStart != leftEnd)
    {
        _leftNode = std::unique_ptr<KDNode<T>>(new KDNode<T>);
        _leftNode->build(data, leftStart, leftEnd, _axis);
    }

    if (rightStart != rightEnd)
    {
        _rightNode = std::unique_ptr<KDNode<T>>(new KDNode<T>);
        _rightNode->build(data, rightStart, rightEnd, _axis);
    }
}


template <class T>
void KDNode<T>::nearestNeighbor(const std::vector<T>& queryPoint,
    IndexedPoint<T>& best,
    double& bestSqrDistance)
{
    using namespace std;

    assert(_location.point.size() == queryPoint.size());

    // If the point at this node is closer than our current best, make
    // it the best
    T distance = squaredDistance<T>(_location.point, queryPoint);
    if (distance < bestSqrDistance)
    {
        bestSqrDistance = distance;
        best = _location;
    }

    // Recursively search the nodes beneath this one
    if (_leftNode != nullptr && queryPoint[_axis] <= _location.point[_axis])
    {
        _leftNode->nearestNeighbor(queryPoint, best, bestSqrDistance);
    }
    else if (_rightNode != nullptr)
    {
        _rightNode->nearestNeighbor(queryPoint, best, bestSqrDistance);
    }

    // If the this node's splitting plane is within the radius of the
    // best distance hypersphere, we need to examine the node on the
    // other side of plane, in case it is actually closer
    const double hypersphereDist = queryPoint[_axis] - _location.point[_axis];
    if (hypersphereDist * hypersphereDist <= bestSqrDistance)
    {
        if (_rightNode != nullptr && queryPoint[_axis] <= _location.point[_axis])
        {
            _rightNode->nearestNeighbor(queryPoint, best, bestSqrDistance);
        }
        else if (_leftNode != nullptr)
        {
            _leftNode->nearestNeighbor(queryPoint, best, bestSqrDistance);
        }
    }
}