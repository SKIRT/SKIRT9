#include "SearchTree.hpp"
#include "FatalError.hpp"

using Node = SearchTree::Node;

namespace
{
    // returns the square of its argument
    inline double sqr(double x)
    {
        return x * x;
    }

    // function to compare two points according to the specified axis (0,1,2)
    inline bool lessthan(Vec p1, Vec p2, int axis)
    {
        switch (axis)
        {
            case 0:  // split on x
                if (p1.x() < p2.x()) return true;
                if (p1.x() > p2.x()) return false;
                if (p1.y() < p2.y()) return true;
                if (p1.y() > p2.y()) return false;
                if (p1.z() < p2.z()) return true;
                return false;
            case 1:  // split on y
                if (p1.y() < p2.y()) return true;
                if (p1.y() > p2.y()) return false;
                if (p1.z() < p2.z()) return true;
                if (p1.z() > p2.z()) return false;
                if (p1.x() < p2.x()) return true;
                return false;
            case 2:  // split on z
                if (p1.z() < p2.z()) return true;
                if (p1.z() > p2.z()) return false;
                if (p1.x() < p2.x()) return true;
                if (p1.x() > p2.x()) return false;
                if (p1.y() < p2.y()) return true;
                return false;
            default:  // this should never happen
                return false;
        }
    }
}

////////////////////////////////////////////////////////////////////

SearchTree::~SearchTree()
{
    delete _root;
}

////////////////////////////////////////////////////////////////////

SearchTree::SearchTree() {}

////////////////////////////////////////////////////////////////////

void SearchTree::buildTree(const vector<Vec>& nodePositions)
{
    vector<int> indices(nodePositions.size());
    std::iota(indices.begin(), indices.end(), 0);  // just use loop for clarity?
    _root = buildTree(nodePositions, indices.begin(), indices.end(), 0);
}

////////////////////////////////////////////////////////////////////

void SearchTree::buildTree(const vector<Vec>& nodePositions, vector<int>& indices)
{
    _root = buildTree(nodePositions, indices.begin(), indices.end(), 0);
}

////////////////////////////////////////////////////////////////////

Node* SearchTree::buildTree(const vector<Vec>& nodePositions, const vector<int>::iterator first,
                            const vector<int>::iterator last, int depth) const
{
    auto length = last - first;
    if (length > 0)
    {
        auto median = length >> 1;
        std::nth_element(first, first + median, last, [&nodePositions, depth](int m1, int m2) {
            return m1 != m2 && lessthan(nodePositions[m1], nodePositions[m2], depth % 3);
        });
        return new Node(*(first + median), depth, buildTree(nodePositions, first, first + median, depth + 1),
                        buildTree(nodePositions, first + median + 1, last, depth + 1));
    }
    return nullptr;
}

////////////////////////////////////////////////////////////////////F

int SearchTree::nearest(Vec& bfr, const vector<Vec>& nodePositions) const
{
    return _root->nearest(bfr, nodePositions)->m();
}

////////////////////////////////////////////////////////////////////F

Node* SearchTree::root() const
{
    return _root;
}

////////////////////////////////////////////////////////////////////

Node::Node(int m, int depth, Node* left, Node* right) : _m(m), _axis(depth % 3), _up(0), _left(left), _right(right)
{
    if (_left) _left->setParent(this);
    if (_right) _right->setParent(this);
}

////////////////////////////////////////////////////////////////////

Node::~Node()
{
    delete _left;
    delete _right;
}

////////////////////////////////////////////////////////////////////

void Node::setParent(Node* up)
{
    _up = up;
}

////////////////////////////////////////////////////////////////////

int Node::m() const
{
    return _m;
}

////////////////////////////////////////////////////////////////////

Node* Node::up() const
{
    return _up;
}

////////////////////////////////////////////////////////////////////

Node* Node::left() const
{
    return _left;
}

////////////////////////////////////////////////////////////////////

Node* Node::right() const
{
    return _right;
}

////////////////////////////////////////////////////////////////////

Node* Node::child(Vec bfr, const vector<Vec>& nodePositions) const
{
    return lessthan(bfr, nodePositions[_m], _axis) ? _left : _right;
}

////////////////////////////////////////////////////////////////////

Node* Node::otherChild(Vec bfr, const vector<Vec>& nodePositions) const
{
    return lessthan(bfr, nodePositions[_m], _axis) ? _right : _left;
}

////////////////////////////////////////////////////////////////////

const Vec& Node::position(const vector<Vec>& nodePositions) const
{
    return nodePositions[_m];
}

////////////////////////////////////////////////////////////////////

double Node::squaredDistanceToSplitPlane(Vec bfr, const vector<Vec>& nodePositions) const
{
    switch (_axis)
    {
        case 0:  // split on x
            return sqr(nodePositions[_m].x() - bfr.x());
        case 1:  // split on y
            return sqr(nodePositions[_m].y() - bfr.y());
        case 2:  // split on z
            return sqr(nodePositions[_m].z() - bfr.z());
        default:  // this should never happen
            return 0;
    }
}

////////////////////////////////////////////////////////////////////

Node* Node::nearest(Vec bfr, const vector<Vec>& nodePositions)
{
    // recursively descend the tree until a leaf node is reached, going left or right depending on
    // whether the specified point is less than or greater than the current node in the split dimension
    Node* current = this;
    Node* child = current->child(bfr, nodePositions);
    while (child)
    {
        current = child;
        child = current->child(bfr, nodePositions);
    }

    // unwind the recursion, looking for the nearest node while climbing up
    Node* best = current;
    Vec bestPos = best->position(nodePositions);
    double bestSD = (bestPos - bfr).norm2();
    while (true)
    {
        // if the current node is closer than the current best, then it becomes the current best
        Vec currentPos = current->position(nodePositions);
        double currentSD = (currentPos - bfr).norm2();
        if (currentSD < bestSD)
        {
            best = current;
            bestSD = currentSD;
        }

        // if there could be points on the other side of the splitting plane for the current node
        // that are closer to the search point than the current best, then ...
        double splitSD = current->squaredDistanceToSplitPlane(bfr, nodePositions);
        if (splitSD < bestSD)
        {
            // move down the other branch of the tree from the current node looking for closer points,
            // following the same recursive process as the entire search
            Node* other = current->otherChild(bfr, nodePositions);
            if (other)
            {
                Node* otherBest = other->nearest(bfr, nodePositions);
                Vec otherBestPos = otherBest->position(nodePositions);
                double otherBestSD = (otherBestPos - bfr).norm2();
                if (otherBestSD < bestSD)
                {
                    best = otherBest;
                    bestSD = otherBestSD;
                }
            }
        }

        // move up to the parent until we meet the top node
        if (current == this) break;
        current = current->up();
    }
    return best;
}

////////////////////////////////////////////////////////////////////