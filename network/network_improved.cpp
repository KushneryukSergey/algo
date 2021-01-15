#include <iostream>
#include <vector>
#include <queue>
#include <list>
#include <string>
#include <algorithm>
#include <tuple>

using std::string;
using std::vector;
using std::queue;
using std::list;


enum EDGE_TYPE : bool {
    DIRECT_EDGE = false,
    REVERSE_EDGE = true
};

// * All vertices have numbers from 0 to N - 1

// * get_target() and get_source() functions return indices of source and target

// * size() function returns number of vertices

// * begin(vertex, DIRECT_EDGE) should return forward iterator over ALL edges
//   which start in vertex

// * begin(vertex, REVERSE_EDGE) should return forward iterator over ALL edges
//   which end in vertex

// * Same is required for end(vertex, DIRECT_EDGE)
template<typename FlowType>
class NetworkInterface {
public:
    [[nodiscard]] virtual const size_t& get_source() const = 0;
    
    [[nodiscard]] virtual const size_t& get_target() const = 0;
    
    [[nodiscard]] virtual size_t size() const = 0;
    
    class EdgeIteratorInterface;
    
    class iterator;
    
    [[nodiscard]] virtual iterator begin(const size_t&, EDGE_TYPE) = 0;
    
    [[nodiscard]] virtual iterator end(const size_t&, EDGE_TYPE) = 0;
};

template<typename FlowType>
class NetworkInterface<FlowType>::iterator {
public:
    iterator() = delete;
    
    explicit iterator(EdgeIteratorInterface*);
    iterator(iterator&&) noexcept ;
    
    bool operator!=(iterator const& other) const;
    
    bool operator==(iterator const& other) const;
    
    ~iterator();
    
    iterator& operator++();
    
    iterator& operator=(iterator&&);
    
    FlowType get_free_capacity() const;
    
    [[nodiscard]] const size_t& get_start() const;
    
    [[nodiscard]] const size_t& get_finish() const;
    
    void push(const FlowType&) const;
    
    void push_max() const;

private:
    EdgeIteratorInterface* _pointer;
};

template<typename FlowType>
class NetworkInterface<FlowType>::EdgeIteratorInterface {
public:
    EdgeIteratorInterface() = delete;
    
    explicit EdgeIteratorInterface(NetworkInterface<FlowType>&);
    
    virtual ~EdgeIteratorInterface() = default;
    
    [[nodiscard]] virtual bool is_equal(const EdgeIteratorInterface& other) const = 0;
    
    virtual EdgeIteratorInterface& operator++() = 0;
    
    [[nodiscard]] virtual FlowType get_free_capacity() const = 0;
    
    [[nodiscard]] virtual const size_t& get_start() const = 0;
    
    [[nodiscard]] virtual const size_t& get_finish() const = 0;
    
    virtual void push(const FlowType&) const = 0;
    
    virtual void push_max() const;

protected:
    NetworkInterface& _network;
};

template<typename FlowType>
NetworkInterface<FlowType>::iterator::iterator(
        NetworkInterface<FlowType>::EdgeIteratorInterface* pointer):
        _pointer(pointer) {}

template<typename FlowType>
bool NetworkInterface<FlowType>::iterator::operator!=(
        const NetworkInterface<FlowType>::iterator& other) const {
    return !(*this == other);
}

template<typename FlowType>
bool NetworkInterface<FlowType>::iterator::operator==(
        const NetworkInterface<FlowType>::iterator& other) const {
    return _pointer->is_equal(*other._pointer);
}

template<typename FlowType>
typename NetworkInterface<FlowType>::iterator&
NetworkInterface<FlowType>::iterator::operator++() {
    ++(*_pointer);
    return *this;
}

template<typename FlowType>
typename NetworkInterface<FlowType>::iterator&
NetworkInterface<FlowType>::iterator::operator=(
        NetworkInterface<FlowType>::iterator&& other) {
    //delete _pointer;
    _pointer = other._pointer;
    other._pointer = nullptr;
    return *this;
}

template<typename FlowType>
FlowType NetworkInterface<FlowType>::iterator::get_free_capacity() const {
    return _pointer->get_free_capacity();
}

template<typename FlowType>
const size_t& NetworkInterface<FlowType>::iterator::get_start() const {
    return _pointer->get_start();
}

template<typename FlowType>
const size_t& NetworkInterface<FlowType>::iterator::get_finish() const {
    return _pointer->get_finish();
}

template<typename FlowType>
void NetworkInterface<FlowType>::iterator::push(const FlowType& value) const {
    return _pointer->push(value);
}

template<typename FlowType>
void NetworkInterface<FlowType>::iterator::push_max() const {
    return _pointer->push_max();
}

template<typename FlowType>
NetworkInterface<FlowType>::iterator::~iterator() {
    delete _pointer;
}

template<typename FlowType>
NetworkInterface<FlowType>::iterator::iterator(
        NetworkInterface<FlowType>::iterator&& other) noexcept:
            _pointer(other._pointer){
    other._pointer = nullptr;
}


template<typename FlowType>
void NetworkInterface<FlowType>::EdgeIteratorInterface::push_max() const {
    push(get_free_capacity());
}

template<typename FlowType>
NetworkInterface<FlowType>::EdgeIteratorInterface::EdgeIteratorInterface(
        NetworkInterface<FlowType>& network): _network(network) {}


template<typename FlowType>
class Network : public NetworkInterface<FlowType> {
private:
    typedef NetworkInterface<FlowType> Base;
    typedef typename Base::iterator iterator;
public:
    // Network own methods
    
    Network() = delete;
    
    explicit Network(const size_t&, const size_t&, const size_t&);
    
    [[nodiscard]] const size_t& get_source() const override;
    
    [[nodiscard]] const size_t& get_target() const override;
    
    [[nodiscard]] size_t size() const override;
    
    
    // Methods connected with Edge
    class Edge;
    
    void add_edge(const size_t&, const size_t&, const FlowType&);
    
    [[nodiscard]] const Edge& get_edge(const size_t&) const;
    
    void push(const size_t&, const FlowType&);
    
    
    // Methods connected with EdgeIterator
    class EdgeIterator;
    
    iterator begin(const size_t&, EDGE_TYPE) override;
    
    iterator end(const size_t&, EDGE_TYPE) override;

private:
    size_t _source;
    size_t _target;
    
    vector<Edge> _edges;
    vector<vector<size_t>> _graph;
    vector<vector<size_t>> _rev_graph;
};


// -----------------------------------------
//              Edge methods
// -----------------------------------------
template<typename FlowType>
class Network<FlowType>::Edge {
public:
    Edge() = delete;
    
    Edge(const size_t&, const size_t&, const FlowType&);
    
    [[nodiscard]] const size_t& get_finish() const;
    
    [[nodiscard]] const size_t& get_start() const;
    
    void push(const FlowType&);
    
    FlowType get_free_capacity() const;

private:
    size_t _start, _finish;
    FlowType _capacity;
    FlowType _flow = FlowType(0);
};

template<typename FlowType>
Network<FlowType>::Edge::Edge(const size_t& start,
                              const size_t& finish,
                              const FlowType& capacity):
        _start(start),
        _finish(finish),
        _capacity(capacity) {}

template<typename FlowType>
void Network<FlowType>::Edge::push(const FlowType& value) {
    _flow += value;
}

template<typename FlowType>
FlowType Network<FlowType>::Edge::get_free_capacity() const {
    return _capacity - _flow;
}

template<typename FlowType>
const size_t& Network<FlowType>::Edge::get_finish() const {
    return _finish;
}

template<typename FlowType>
const size_t& Network<FlowType>::Edge::get_start() const {
    return _start;
}




// -----------------------------------------
//           EdgeIterator methods
// -----------------------------------------

template<typename FlowType>
class Network<FlowType>::EdgeIterator :
        public NetworkInterface<FlowType>::EdgeIteratorInterface {
private:
    typedef typename NetworkInterface<FlowType>::EdgeIteratorInterface Base;
    vector<size_t>::iterator it;
public:
    EdgeIterator() = delete;
    
    explicit EdgeIterator(const size_t&, Network<FlowType>&, const bool&);
    
    explicit EdgeIterator(const vector<size_t>::iterator&, Network<FlowType>&);
    
    ~EdgeIterator() override = default;
    
    [[nodiscard]] bool is_equal(const Base& other) const override;
    
    EdgeIterator& operator++() override;
    
    EdgeIterator& operator=(const EdgeIterator&);
    
    [[nodiscard]] FlowType get_free_capacity() const override;
    
    [[nodiscard]] const size_t& get_start() const override;
    
    [[nodiscard]] const size_t& get_finish() const override;
    
    void push(const FlowType&) const override;
};

template<typename FlowType>
Network<FlowType>::EdgeIterator::EdgeIterator(
        const size_t& vertex_index,
        Network<FlowType>& network,
        const bool& reverse):
        NetworkInterface<FlowType>::EdgeIteratorInterface(network),
        it((reverse ? network._rev_graph : network._graph)[vertex_index].
                                                                   begin()) {}

template<typename FlowType>
Network<FlowType>::EdgeIterator::EdgeIterator(
        const vector<size_t>::iterator& it,
        Network& network):
        NetworkInterface<FlowType>::EdgeIteratorInterface(network),
        it(it) {}


template<typename FlowType>
typename Network<FlowType>::EdgeIterator&
Network<FlowType>::EdgeIterator::operator++() {
    ++it;
    return *this;
}

template<typename FlowType>
typename Network<FlowType>::EdgeIterator&
Network<FlowType>::EdgeIterator::operator=(
        const EdgeIterator& other) {
    it = other.it;
    return *this;
}

template<typename FlowType>
FlowType Network<FlowType>::EdgeIterator::get_free_capacity() const {
    return dynamic_cast<Network<FlowType>&>(Base::_network).
                                get_edge(*it).get_free_capacity();
}

template<typename FlowType>
const size_t& Network<FlowType>::EdgeIterator::get_start() const {
    return dynamic_cast<Network<FlowType>&>(Base::_network).
                                        get_edge(*it).get_start();
}

template<typename FlowType>
const size_t& Network<FlowType>::EdgeIterator::get_finish() const {
    return dynamic_cast<Network<FlowType>&>(Base::_network).
                                        get_edge(*it).get_finish();
}

template<typename FlowType>
void Network<FlowType>::EdgeIterator::push(const FlowType& value) const {
    dynamic_cast<Network<FlowType>&>(Base::_network).push(*it, value);
}

template<typename FlowType>
bool Network<FlowType>::EdgeIterator::is_equal(const Base& other) const {
    return
    it == dynamic_cast<const Network<FlowType>::EdgeIterator&>(other).it;
}





// -----------------------------------------
//            Network Methods
// -----------------------------------------

template<typename FlowType>
Network<FlowType>::Network(const size_t& source,
                           const size_t& target,
                           const size_t& vertex_number):
        _source(source),
        _target(target) {
    _graph.resize(vertex_number);
    _rev_graph.resize(vertex_number);
}

template<typename FlowType>
const size_t& Network<FlowType>::get_source() const {
    return _source;
}

template<typename FlowType>
const size_t& Network<FlowType>::get_target() const {
    return _target;
}

template<typename FlowType>
size_t Network<FlowType>::size() const {
    return _graph.size();
}


template<typename FlowType>
void Network<FlowType>::add_edge(const size_t& start,
                                 const size_t& finish,
                                 const FlowType& capacity) {
    _edges.emplace_back(start, finish, capacity);
    _edges.emplace_back(finish, start, 0);
    
    _graph[start].push_back(_edges.size() - 2);
    _graph[finish].push_back(_edges.size() - 1);
    
    _rev_graph[finish].push_back(_edges.size() - 2);
    _rev_graph[start].push_back(_edges.size() - 1);
}

template<typename FlowType>
const typename Network<FlowType>::Edge& Network<FlowType>::get_edge(
        const size_t& edge_index) const {
    return _edges[edge_index];
}

template<typename FlowType>
void
Network<FlowType>::push(const size_t& edge_index, const FlowType& value) {
    _edges[edge_index].push(value);
    _edges[edge_index ^ 1u].push(-value);
}


template<typename FlowType>
typename Network<FlowType>::iterator Network<FlowType>::begin(
        const size_t& vertex_index,
        EDGE_TYPE type) {
    return iterator(new EdgeIterator(vertex_index, *this, type));
}

template<typename FlowType>
typename Network<FlowType>::iterator Network<FlowType>::end(
        const size_t& vertex_index,
        EDGE_TYPE type) {
    auto pointer = new EdgeIterator(
            (type ? _rev_graph : _graph)[vertex_index].end(),
            *this);
    return iterator(pointer);
}














// -----------------------------------------
//         Network Max Flow Finders
// -----------------------------------------
/* *
 * Class for solving classical "Finding maximal
 * flow in an arbitrary network" problem
 */

template<typename FlowType>
class NetworkMaxFlowFinder {
protected:
    NetworkInterface<FlowType>& network;
    size_t source;
    size_t target;
public:
    explicit NetworkMaxFlowFinder(NetworkInterface<FlowType>&);
    
    //--------------------------------------
    // Your majesty, Find-Max-Flow Function
    //--------------------------------------
    virtual FlowType find_max_flow() = 0;
};

template<typename FlowType>
NetworkMaxFlowFinder<FlowType>::NetworkMaxFlowFinder(
        NetworkInterface<FlowType>& network):
        network(network),
        source(network.get_source()),
        target(network.get_target()) {}




// -----------------------------------------------------
//      "Blocking Flow" Technology Max Flow Finder
// -----------------------------------------------------
/* *
 * Parent class for all max-flow finders using blocking
 * flow technology. It's necessary to override
 * find_blocking_flow and set_up functions:
 * _set_up_before(_after) - function, that implements
 * redefining or finding additional information for
 * algorithm you are using before and after refinding
 * levels (by default they do nothing)
 *
 * ATTENTION:
 * You should update answer by "_update_answer" function
 * in "find_blocking_flow", because this parent class
 * return "answer" value in "find_max_flow" function
 *
 *
 * This parent class use layers/levels technology, which
 * means, it find the smallest distances from source to
 * every vertex and find them before every iteration of
 * searching of blocking flow and save it in "_level"
 * array
 */

template<typename FlowType>
class BlockingFlowMaxFlowFinder : protected NetworkMaxFlowFinder<FlowType> {
private:
    FlowType _answer = FlowType(0);
protected:
    vector<size_t> _level;
    
    void _update_answer(const FlowType&);
    
    void _calculate_level();
    
    void _prepare_for_new_iter();
    
    virtual void _find_blocking_flow() = 0;
    
    virtual void _set_up_before() {};
    
    virtual void _set_up_after() {};

public:
    explicit BlockingFlowMaxFlowFinder(NetworkInterface<FlowType>&);
    
    //--------------------------------------
    // Your majesty, Find-Max-Flow Function
    //--------------------------------------
    FlowType find_max_flow() override;
    
};

template<typename FlowType>
BlockingFlowMaxFlowFinder<FlowType>::BlockingFlowMaxFlowFinder(
        NetworkInterface<FlowType>& network):
        NetworkMaxFlowFinder<FlowType>(network) {
    _level.resize(network.size(), SIZE_MAX);
}


template<typename FlowType>
FlowType BlockingFlowMaxFlowFinder<FlowType>::find_max_flow() {
    this->_prepare_for_new_iter();
    while (_level[this->target] != SIZE_MAX) {
        this->_find_blocking_flow();
        this->_prepare_for_new_iter();
    }
    return _answer;
}


template<typename FlowType>
void BlockingFlowMaxFlowFinder<FlowType>::_calculate_level() {
    vector<bool> visited(this->network.size());
    queue<size_t> to_visit;
    to_visit.push(this->source);
    _level[this->source] = 0;
    
    while (!to_visit.empty()) {
        size_t vertex = to_visit.front();
        to_visit.pop();
        if (!visited[vertex]) {
            visited[vertex] = true;
            for (auto vertex_edge_it = this->network.begin(vertex,
                                                           DIRECT_EDGE);
                 vertex_edge_it != this->network.end(vertex,
                                                     DIRECT_EDGE);
                 ++vertex_edge_it) {
                if (vertex_edge_it.get_free_capacity() != FlowType(0)) {
                    size_t finish = vertex_edge_it.get_finish();
                    if (!visited[finish]) {
                        _level[finish] = _level[vertex] + 1;
                        to_visit.push(finish);
                    }
                }
            }
        }
    }
}

template<typename FlowType>
void BlockingFlowMaxFlowFinder<FlowType>::_prepare_for_new_iter() {
    _set_up_before();
    _level.assign(_level.size(), SIZE_MAX);
    _calculate_level();
    _set_up_after();
}

template<typename FlowType>
void BlockingFlowMaxFlowFinder<FlowType>::_update_answer(
        const FlowType& value) {
    _answer += value;
}


// -----------------------------------------------------
//   Malhotra-Kumar-Maheshwari Network Max Flow Finder
// -----------------------------------------------------
/* *
 * Find maximum flow in arbitary network for O(V^3)
 * Using concept of blocking flows, it finds blocking flow
 * in O(V^2) and make O(V) iterations of searching blocking
 * flow (O(V^3) as a whole)
 */

template<typename FlowType>
class MKMMaxFlowFinder : public BlockingFlowMaxFlowFinder<FlowType> {
    typedef NetworkMaxFlowFinder<FlowType> Finder; // shortcut for SearchEngine to get
    // shorter expression to access network field
private:
    vector<FlowType> _potential_in;
    vector<FlowType> _potential_out;
    vector<bool> _deleted;
    vector<typename NetworkInterface<FlowType>::iterator> _direct_iterator,
            _reverse_iterator;
    
    
    void _calculate_potential();
    
    FlowType _get_potential(const size_t&);
    
    size_t _find_min_potential_vertex();
    
    void _find_blocking_flow() override;
    
    void _set_up_after() override;
    
    void _push_to_target(const size_t&, const FlowType&);
    
    void _push_to_source(const size_t&, const FlowType&);
    
    bool _is_correct(const typename NetworkInterface<FlowType>::iterator&);
    
    void _delete_vertex(const size_t&);
    
    void _change_potential_by_edge(
            const typename NetworkInterface<FlowType>::iterator&,
            const bool& push_max = true,
            const FlowType& value = FlowType(0));

public:
    MKMMaxFlowFinder(const MKMMaxFlowFinder&) = delete;
    
    MKMMaxFlowFinder(MKMMaxFlowFinder&&) = delete;
    
    MKMMaxFlowFinder& operator=(const MKMMaxFlowFinder&) = delete;
    
    MKMMaxFlowFinder& operator=(MKMMaxFlowFinder&&) = delete;
    
    explicit MKMMaxFlowFinder(NetworkInterface<FlowType>&);
};

template<typename FlowType>
MKMMaxFlowFinder<FlowType>::MKMMaxFlowFinder(
        NetworkInterface<FlowType>& network):
        BlockingFlowMaxFlowFinder<FlowType>(network) {
    _potential_in.resize(network.size(), FlowType(0));
    _potential_out.resize(network.size(), FlowType(0));
    _deleted.resize(network.size(), false);
    for (size_t i = 0; i < this->network.size(); ++i) {
        _direct_iterator.push_back(
                this->network.begin(i, DIRECT_EDGE));
        _reverse_iterator.push_back(
                this->network.begin(i, REVERSE_EDGE));
    }
}

template<typename FlowType>
void MKMMaxFlowFinder<FlowType>::_calculate_potential() {
    for (size_t vertex = 0; vertex < this->network.size(); ++vertex) {
        for (auto vertex_edge_it = this->network.begin(
                vertex, DIRECT_EDGE);
             vertex_edge_it != this->network.end(
                     vertex, DIRECT_EDGE);
             ++vertex_edge_it) {
            if (_is_correct(vertex_edge_it)) {
                FlowType capacity = vertex_edge_it.get_free_capacity();
                _potential_out[vertex] += capacity;
                _potential_in[vertex_edge_it.get_finish()] += capacity;
            }
        }
    }
}


template<typename FlowType>
FlowType MKMMaxFlowFinder<FlowType>::_get_potential(const size_t& vertex) {
    if (vertex == this->source) {
        return _potential_out[vertex];
    }
    if (vertex == this->target) {
        return _potential_in[vertex];
    }
    return std::min(_potential_in[vertex], _potential_out[vertex]);
}

template<typename FlowType>
size_t MKMMaxFlowFinder<FlowType>::_find_min_potential_vertex() {
    size_t min_potential_vertex = this->source;
    for (size_t i = 0; i < this->network.size(); ++i) {
        if (!_deleted[i] &&
            _get_potential(i) < _get_potential(min_potential_vertex)) {
            min_potential_vertex = i;
        }
    }
    return min_potential_vertex;
}

template<typename FlowType>
void MKMMaxFlowFinder<FlowType>::_find_blocking_flow() {
    while (!_deleted[this->source]) {
        size_t min_potential_vertex = _find_min_potential_vertex();
        FlowType min_potential = _get_potential(min_potential_vertex);
        if (min_potential == FlowType(0)) {
            _delete_vertex(min_potential_vertex);
        } else {
            this->_update_answer(_get_potential(min_potential_vertex));
            _push_to_source(min_potential_vertex, min_potential);
            _push_to_target(min_potential_vertex, min_potential);
        }
    }
}

template<typename FlowType>
void MKMMaxFlowFinder<FlowType>::_set_up_after() {
    _deleted.assign(_deleted.size(), false);
    _potential_in.assign(_potential_in.size(), 0);
    _potential_out.assign(_potential_out.size(), 0);
    for (size_t i = 0; i < this->network.size(); ++i) {
        _direct_iterator[i] = this->network.begin(i,
                                                  DIRECT_EDGE);
        _reverse_iterator[i] = this->network.begin(i,
                                                   REVERSE_EDGE);
    }
    _calculate_potential();
}

template<typename FlowType>
void MKMMaxFlowFinder<FlowType>::_push_to_target(
        const size_t& min_potential_vertex,
        const FlowType& potential) {
    vector<FlowType> to_be_pushed(this->network.size());
    queue<size_t> blocking_flow_vertexes;
    blocking_flow_vertexes.push(min_potential_vertex);
    to_be_pushed[min_potential_vertex] = potential;
    
    
    while (!blocking_flow_vertexes.empty()) {
        size_t vertex = blocking_flow_vertexes.front();
        blocking_flow_vertexes.pop();
        while (to_be_pushed[vertex] != FlowType(0) && vertex != this->target) {
            auto& edge_it = _direct_iterator[vertex];
            if (_is_correct(edge_it)) {
                if (to_be_pushed[edge_it.get_finish()] == FlowType(0)) {
                    blocking_flow_vertexes.push(edge_it.get_finish());
                }
                if (edge_it.get_free_capacity() > to_be_pushed[vertex]) {
                    _change_potential_by_edge(edge_it,
                                              false,
                                              to_be_pushed[vertex]);
                    to_be_pushed[edge_it.get_finish()] += to_be_pushed[vertex];
                    edge_it.push(to_be_pushed[vertex]);
                    to_be_pushed[vertex] = 0;
                    break;
                }
                _change_potential_by_edge(edge_it);
                to_be_pushed[edge_it.get_finish()] += edge_it.get_free_capacity();
                to_be_pushed[vertex] -= edge_it.get_free_capacity();
                edge_it.push_max();
            }
            ++_direct_iterator[vertex];
        }
    }
}

template<typename FlowType>
void MKMMaxFlowFinder<FlowType>::_push_to_source(const size_t& min_potential_vertex,
                                                 const FlowType& potential) {
    vector<FlowType> to_be_pushed(this->network.size());
    queue<size_t> blocking_flow_vertexes;
    blocking_flow_vertexes.push(min_potential_vertex);
    to_be_pushed[min_potential_vertex] = potential;
    
    
    while (!blocking_flow_vertexes.empty()) {
        size_t vertex = blocking_flow_vertexes.front();
        blocking_flow_vertexes.pop();
        while (to_be_pushed[vertex] != FlowType(0) && vertex != this->source) {
            auto& edge_it = _reverse_iterator[vertex];
            if (_is_correct(edge_it)) {
                if (to_be_pushed[edge_it.get_start()] == FlowType(0)) {
                    blocking_flow_vertexes.push(edge_it.get_start());
                }
                if (edge_it.get_free_capacity() > to_be_pushed[vertex]) {
                    _change_potential_by_edge(edge_it,
                                              false,
                                              to_be_pushed[vertex]);
                    to_be_pushed[edge_it.get_start()] += to_be_pushed[vertex];
                    edge_it.push(to_be_pushed[vertex]);
                    to_be_pushed[vertex] = 0;
                    break;
                }
                _change_potential_by_edge(edge_it);
                to_be_pushed[edge_it.get_start()] += edge_it.get_free_capacity();
                to_be_pushed[vertex] -= edge_it.get_free_capacity();
                edge_it.push_max();
            }
            ++_reverse_iterator[vertex];
        }
    }
}

template<typename FlowType>
bool MKMMaxFlowFinder<FlowType>::_is_correct(
        const typename NetworkInterface<FlowType>::iterator& it) {
    return this->_level[it.get_start()] + 1 == this->_level[it.get_finish()] &&
           !_deleted[it.get_start()] &&
           !_deleted[it.get_finish()] &&
           it.get_free_capacity() > FlowType(0);
}


template<typename FlowType>
void MKMMaxFlowFinder<FlowType>::_delete_vertex(const size_t& vertex) { // delete if potential is 0
    if (_get_potential(vertex) != FlowType(0)) {
        return;
    }
    auto iterator_type = _potential_in[vertex] == FlowType(0) ?
                         DIRECT_EDGE :
                         REVERSE_EDGE;
    
    
    for (auto vertex_edge_it = this->network.begin(vertex, iterator_type);
         vertex_edge_it != this->network.end(vertex, iterator_type);
         ++vertex_edge_it) {
        if (_is_correct(vertex_edge_it)) {
            _change_potential_by_edge(vertex_edge_it);
        }
    }
    _deleted[vertex] = true;
}

template<typename FlowType>
void MKMMaxFlowFinder<FlowType>::_change_potential_by_edge(
        const typename NetworkInterface<FlowType>::iterator& it,
        const bool& push_max,
        const FlowType& value) {
    if (push_max) {
        _potential_in[it.get_finish()] -= it.get_free_capacity();
        _potential_out[it.get_start()] -= it.get_free_capacity();
    } else {
        _potential_in[it.get_finish()] -= value;
        _potential_out[it.get_start()] -= value;
    }
}







// -----------------------------------------------------
//    Goldberg [push-relabel] Network Max Flow Finder
// -----------------------------------------------------
/* *
 * Alternative Goldberg technology of searching maximum flow
 * in a network by transforming pre-flow into max flow
 * Point-blank algorithm of "making push while it's possible"
 * works in O(EV^2). It has several improved versions working
 * in O(V^3) using discharge: relabel-to-front, FIFO and Highest
 * label selection rules
 *
 * Push, relabel functions are implemented.
 * In overrided "find_max_flow" function you just should call
 * discharge in correct order
 *
 * ATTENTION: Discharge function added in derived Goldberg
 * Discharge class
 */

template<typename FlowType>
class GoldbergMaxFlowFinder : protected NetworkMaxFlowFinder<FlowType> {
public:
    explicit GoldbergMaxFlowFinder(NetworkInterface<FlowType>&);
    
    //--------------------------------------
    // Your majesty, Find-Max-Flow Function
    //--------------------------------------
    FlowType find_max_flow() override = 0;

protected:
    vector<size_t> _height;
    vector<FlowType> _excess;
    
    bool _is_possible_to_push(
            const typename NetworkInterface<FlowType>::iterator&);
    
    void _push(const typename NetworkInterface<FlowType>::iterator&);
    
    void _relabel(const size_t&);
    
    void _push_from_source();
};

template<typename FlowType>
GoldbergMaxFlowFinder<FlowType>::GoldbergMaxFlowFinder(
        NetworkInterface<FlowType>& network):
        NetworkMaxFlowFinder<FlowType>(network) {
    _height.resize(network.size(), 0);
    _height[network.get_source()] = network.size();
    _excess.resize(network.size(), FlowType(0));
}


template<typename FlowType>
void GoldbergMaxFlowFinder<FlowType>::_push(
        const typename NetworkInterface<FlowType>::iterator& it) {
    FlowType value = std::min(it.get_free_capacity(), _excess[it.get_start()]);
    it.push(value);
    _excess[it.get_finish()] += value;
    _excess[it.get_start()] -= value;
}

template<typename FlowType>
void GoldbergMaxFlowFinder<FlowType>::_relabel(const size_t& vertex) {
    size_t new_height = 2 * this->network.size() + 1;
    for (auto it = this->network.begin(vertex, DIRECT_EDGE);
         it != this->network.end(vertex, DIRECT_EDGE);
         ++it) {
        if (it.get_free_capacity() > FlowType(0)) {
            new_height = std::min(new_height, _height[it.get_finish()] + 1);
        }
    }
    _height[vertex] = new_height;
}

template<typename FlowType>
bool GoldbergMaxFlowFinder<FlowType>::_is_possible_to_push(
        const typename NetworkInterface<FlowType>::iterator& it) {
    return (_height[it.get_finish()] + 1 == _height[it.get_start()]) &&
           (it.get_free_capacity() > FlowType(0));
}

template<typename FlowType>
void GoldbergMaxFlowFinder<FlowType>::_push_from_source() {
    size_t source = this->source;
    auto source_excess = FlowType(0);
    for (auto it = this->network.begin(source, DIRECT_EDGE);
         it != this->network.end(source, DIRECT_EDGE);
         ++it) {
        source_excess += it.get_free_capacity();
        _excess[source] += it.get_free_capacity();
        _push(it);
    }
    _excess[source] -= source_excess;
}


// ---------------------------------------------------------
//  Goldberg Discharge Network Max Flow Finder
// ---------------------------------------------------------
/* *
 * Improved version of Goldberg push-relabel technology.
 * Discharge function added
 */


template<typename FlowType>
class GoldbergDischargeMaxFlowFinder :
        protected GoldbergMaxFlowFinder<FlowType> {
public:
    explicit GoldbergDischargeMaxFlowFinder(NetworkInterface<FlowType>&);
    
    //--------------------------------------
    // Your majesty, Find-Max-Flow Function
    //--------------------------------------
    FlowType find_max_flow() override = 0;

protected:
    bool _discharge(const size_t&);
};

template<typename FlowType>
GoldbergDischargeMaxFlowFinder<FlowType>::GoldbergDischargeMaxFlowFinder(
        NetworkInterface<FlowType>& network):
        GoldbergMaxFlowFinder<FlowType>(network) {}


template<typename FlowType>
bool
GoldbergDischargeMaxFlowFinder<FlowType>::_discharge(const size_t& vertex) {
    bool discharged = false;
    for (auto it = this->network.begin(vertex, DIRECT_EDGE);
         it != this->network.end(vertex, DIRECT_EDGE);
         ++it) {
        if (this->_is_possible_to_push(it)) {
            this->_push(it);
            if (this->_excess[vertex] <= FlowType(0)) {
                discharged = true;
                break;
            }
        }
    }
    return discharged;
}

// ---------------------------------------------------------
//  Goldberg RTF [Relabel-To-Front] Network Max Flow Finder
// ---------------------------------------------------------
/* *
 * Improved version of Goldberg discharge technology,
 * works in O(V^3)
 */

template<typename FlowType>
class GoldbergRTFMaxFlowFinder :
        public GoldbergDischargeMaxFlowFinder<FlowType> {
public:
    explicit GoldbergRTFMaxFlowFinder(NetworkInterface<FlowType>&);
    
    //--------------------------------------
    // Your majesty, Find-Max-Flow Function
    //--------------------------------------
    FlowType find_max_flow() override;
};

template<typename FlowType>
GoldbergRTFMaxFlowFinder<FlowType>::GoldbergRTFMaxFlowFinder(
        NetworkInterface<FlowType>& network):
        GoldbergDischargeMaxFlowFinder<FlowType>(network) {}


template<typename FlowType>
FlowType GoldbergRTFMaxFlowFinder<FlowType>::find_max_flow() {
    this->_push_from_source();
    
    list<size_t> vertexes;
    for (size_t i = 0; i < this->network.size(); ++i) {
        if (i != this->source && i != this->network.get_target()) {
            vertexes.push_back(i);
        }
    }
    
    auto it = vertexes.begin();
    while (it != vertexes.end()) {
        if ((this->_excess[*it] == FlowType(0)) || this->_discharge(*it)) {
            ++it;
        } else {
            this->_relabel(*it);
            size_t moved = *it;
            vertexes.erase(it);
            vertexes.push_front(moved);
            it = vertexes.begin();
        }
    }
    
    return -(this->_excess[this->source]);
}


using std::cout;
using std::cin;
using std::endl;

const int BIG_NUMBER = 1000000;

enum ERROR_HANDLER : bool {
    CORRECT = true,
    INCORRECT = false
};

class MatanSolver {
    Network<int> _network1, _network2;
    size_t _N;
    int _answer;
    int _sum_profit;
    
    void _input_vertices(std::istream& input);
    
    void _input_dependencies(std::istream& input);
    
    template<class Finder1, class Finder2>
    std::pair<int, ERROR_HANDLER> _calculate();

public:
    explicit MatanSolver(const size_t&);
    
    template<class Finder1, class Finder2>
    std::pair<int, ERROR_HANDLER> solve(std::istream& input);
};

void MatanSolver::_input_vertices(std::istream& input) {
    int profit;
    for (size_t i = 0; i < _N; ++i) {
        cin >> profit;
        if (profit >= 0) {
            _network1.add_edge(0, i + 1, profit);
            _network2.add_edge(0, i + 1, profit);
            _sum_profit += profit;
        } else {
            _network1.add_edge(i + 1, _N + 1, -profit);
            _network2.add_edge(i + 1, _N + 1, -profit);
        }
    }
}


void MatanSolver::_input_dependencies(std::istream& input) {
    size_t dependencies;
    size_t current_dependence;
    for (size_t i = 0; i < _N; ++i) {
        cin >> dependencies;
        for (size_t j = 0; j < dependencies; ++j) {
            cin >> current_dependence;
            _network1.add_edge(i + 1, current_dependence, BIG_NUMBER);
            _network2.add_edge(i + 1, current_dependence, BIG_NUMBER);
        }
    }
}

template<class Finder1, class Finder2>
std::pair<int, ERROR_HANDLER> MatanSolver::solve(std::istream& input) {
    _sum_profit = 0;
    _input_vertices(input);
    _input_dependencies(input);
    
    return _calculate<Finder1, Finder2>();
}

template<class Finder1, class Finder2>
std::pair<int, ERROR_HANDLER> MatanSolver::_calculate() {
    Finder1 finder1(_network1);
    Finder2 finder2(_network2);
    
    auto pre_ans1 = finder1.find_max_flow();
    auto pre_ans2 = finder2.find_max_flow();
    
    if (pre_ans1 == pre_ans2)
        return {_sum_profit - pre_ans2, CORRECT};
    else
        return {0, INCORRECT};
}

MatanSolver::MatanSolver(const size_t& N) :
        _network1(0, N + 1, N + 2),
        _network2(0, N + 1, N + 2),
        _N(N) {}

int main() {
#ifdef LOCAL
    freopen("input.txt", "rt", stdin);
    //freopen("output.txt", "wt", stdout);
#endif
    
    size_t N;
    cin >> N;
    MatanSolver solver(N);
    int answer;
    ERROR_HANDLER error_check;
    std::tie(answer, error_check) = solver.solve<MKMMaxFlowFinder<int>,
            GoldbergRTFMaxFlowFinder<int>>(cin);
    if (!error_check)
        return 1;
    cout << answer;
    return 0;
}
