#include <iostream>
#include <vector>
#include <queue>
#include <list>
#include <string>
#include <algorithm>

using std::string;
using std::vector;
using std::queue;
using std::list;


// virtual classes
enum EDGE_TYPE : bool {
    DIRECT_EDGE = false,
    REVERSE_EDGE = true
};


template<typename FlowType>
class Network {
public:
    // Network own methods
    
    Network() = delete;
    explicit Network(const size_t&, const size_t&, const size_t&);
    
    [[nodiscard]] const size_t& get_source() const;
    [[nodiscard]] const size_t& get_target() const;
    size_t size();
    
    // Methods connected with Edge
    class Edge;
    void add_edge(const size_t&, const size_t&, const FlowType&);
    const Edge& get_edge(const size_t&) const;
    void push(const size_t&, const FlowType&);
    
    // Methods connected with EdgeIterator
    class EdgeIterator;
    EdgeIterator begin(const size_t&, EDGE_TYPE);
    EdgeIterator end(const size_t&, EDGE_TYPE);
    
    void find_reachable(vector<bool>&);
    
    void push_to_source(const size_t& v, const FlowType& value);
    void push_to_target(const size_t& v, const FlowType& value);
    
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
        public std::iterator<std::forward_iterator_tag,
                typename Network<FlowType>::Edge> {
public:
    EdgeIterator() = delete;
    explicit EdgeIterator(const size_t&, Network<FlowType>&, const bool&);
    explicit EdgeIterator(const vector<size_t>::iterator&, Network<FlowType>&);
    
    bool operator!=(EdgeIterator const& other) const;
    bool operator==(EdgeIterator const& other) const;
    
    EdgeIterator& operator++();
    EdgeIterator& operator=(const EdgeIterator&);
    
    FlowType get_free_capacity() const;
    FlowType get_flow() const;
    [[nodiscard]] const size_t& get_start() const;
    [[nodiscard]] const size_t& get_finish() const;
    void push(const FlowType&) const;
    void push_max() const;

    vector<size_t>::iterator it;
    Network& network;
};

template<typename FlowType>
Network<FlowType>::EdgeIterator::EdgeIterator(const size_t& vertex_index,
                                              Network<FlowType>& network,
                                              const bool& reverse):
        it((reverse ? network._rev_graph : network._graph)[vertex_index].begin()),
        network(network) {}

template<typename FlowType>
Network<FlowType>::EdgeIterator::EdgeIterator(
        const vector<size_t>::iterator& it,
        Network& network):
        it(it), network(network) {}

template<typename FlowType>
bool Network<FlowType>::EdgeIterator::operator==(const EdgeIterator& other) const {
    return it == other.it;
}

template<typename FlowType>
bool Network<FlowType>::EdgeIterator::operator!=(const EdgeIterator& other) const {
    return it != other.it;
}

template<typename FlowType>
typename Network<FlowType>::EdgeIterator& Network<FlowType>::EdgeIterator::operator++() {
    ++it;
    return *this;
}

template<typename FlowType>
typename Network<FlowType>::EdgeIterator& Network<FlowType>::EdgeIterator::operator=(
        const EdgeIterator& other) {
    it = other.it;
    return *this;
}

template<typename FlowType>
FlowType Network<FlowType>::EdgeIterator::get_free_capacity() const {
    return network.get_edge(*it).get_free_capacity();
}

template<typename FlowType>
FlowType Network<FlowType>::EdgeIterator::get_flow() const {
    return network.get_edge(*it)._flow;
}

template<typename FlowType>
const size_t& Network<FlowType>::EdgeIterator::get_start() const {
    return network.get_edge(*it).get_start();
}

template<typename FlowType>
const size_t& Network<FlowType>::EdgeIterator::get_finish() const {
    return network.get_edge(*it).get_finish();
}

template<typename FlowType>
void Network<FlowType>::EdgeIterator::push(const FlowType& value) const {
    network.push(*it, value);
}

template<typename FlowType>
void Network<FlowType>::EdgeIterator::push_max() const {
    push(get_free_capacity());
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
size_t Network<FlowType>::size() {
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
void Network<FlowType>::push(const size_t& edge_index, const FlowType& value) {
    _edges[edge_index].push(value);
    _edges[edge_index ^ 1u].push(-value);
}


template<typename FlowType>
typename Network<FlowType>::EdgeIterator Network<FlowType>::begin(
        const size_t& vertex_index,
        EDGE_TYPE type) {
    return EdgeIterator(vertex_index, *this, type);
}

template<typename FlowType>
typename Network<FlowType>::EdgeIterator Network<FlowType>::end(
        const size_t& vertex_index,
        EDGE_TYPE type) {
    return EdgeIterator((type ? _rev_graph : _graph)[vertex_index].end(),
                        *this);
}

template<typename FlowType>
void Network<FlowType>::find_reachable(vector<bool>& answer) {
    answer.resize(size());
    answer.assign(size(), false);
    
    queue<size_t> q;
    
    q.push(_source);
    answer[_source] = true;
    while(!q.empty()){
        size_t v = q.front();
        q.pop();
        for(auto it = begin(v, DIRECT_EDGE); it != end(v, DIRECT_EDGE); ++it){
            if(!answer[it.get_finish()] && it.get_free_capacity() > FlowType(0)){
                answer[it.get_finish()] = true;
                q.push(it.get_finish());
            }
        }
    }
}

template<typename FlowType>
void Network<FlowType>::push_to_target(const size_t& v, const FlowType& value) {
    if(v == _target) {
        return;
    }
    vector<int> edge(size(), -1);
    
    queue<size_t> q;
    q.push(v);
    edge[v] = -2;
    bool done = false;
    while (!q.empty()) {
        size_t u = q.front();
        q.pop();
        for (auto it = begin(u, DIRECT_EDGE); it != end(u, DIRECT_EDGE); ++it) {
            if (edge[it.get_finish()] == -1 && it.get_flow() >= value) {
                edge[it.get_finish()] = *(it.it);
                q.push(it.get_finish());
                if (it.get_finish() == _target) {
                    done = true;
                    break;
                }
            }
        }
        if (done)
            break;
    }
    int index = edge[_target];
    while(index != -2) {
        push(index, -value);
        index = edge[_edges[index]._start];
    }
}

template<typename FlowType>
void Network<FlowType>::push_to_source(const size_t& v, const FlowType& value) {
    if (v == _source) {
        return;
    }
    vector<int> edge(size(), -1);
    
    queue<size_t> q;
    q.push(v);
    edge[v] = -2;
    bool done = false;
    while (!q.empty()) {
        size_t u = q.front();
        q.pop();
        for (auto it = begin(u, REVERSE_EDGE); it != end(u, REVERSE_EDGE); ++it) {
            if (edge[it.get_start()] == -1 && it.get_flow() >= value) {
                edge[it.get_start()] = *(it.it);
                q.push(it.get_start());
                if (it.get_start() == _source) {
                    done = true;
                    break;
                }
            }
        }
        if (done)
            break;
    }
    int index = edge[_source];
    while (index != -2) {
        push(index, -value);
        index = edge[_edges[index]._finish];
    }
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
    Network<FlowType>& network;
    size_t source;
    size_t target;
public:
    explicit NetworkMaxFlowFinder(Network<FlowType>&);
    
    //--------------------------------------
    // Your majesty, Find-Max-Flow Function
    //--------------------------------------
    virtual FlowType find_max_flow() = 0;
};

template<typename FlowType>
NetworkMaxFlowFinder<FlowType>::NetworkMaxFlowFinder(
        Network<FlowType>& network):
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
    void _precalc_ans();
    virtual void _find_blocking_flow() = 0;
    virtual void _set_up_before() {};
    virtual void _set_up_after() {};
    
public:
    explicit BlockingFlowMaxFlowFinder(Network<FlowType>&);
    
    //--------------------------------------
    // Your majesty, Find-Max-Flow Function
    //--------------------------------------
    FlowType find_max_flow() override;
    
};

template<typename FlowType>
BlockingFlowMaxFlowFinder<FlowType>::BlockingFlowMaxFlowFinder(
        Network<FlowType>& network):
        NetworkMaxFlowFinder<FlowType>(network) {
    _level.resize(network.size(), SIZE_MAX);
}


template<typename FlowType>
FlowType BlockingFlowMaxFlowFinder<FlowType>::find_max_flow() {
    _precalc_ans();
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

template<typename FlowType>
void BlockingFlowMaxFlowFinder<FlowType>::_precalc_ans() {
    _answer = FlowType(0);
    for(auto it = this->network.begin(this->source, DIRECT_EDGE); it != this->network.end(this->source, DIRECT_EDGE); ++it) {
        _answer += it.get_flow();
    }
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
    vector<typename Network<FlowType>::EdgeIterator> _direct_iterator,
            _reverse_iterator;
    
    
    void _calculate_potential();
    FlowType _get_potential(const size_t&);
    size_t _find_min_potential_vertex();
    void _find_blocking_flow() override;
    void _set_up_after() override;
    void _push_to_target(const size_t&, const FlowType&);
    void _push_to_source(const size_t&, const FlowType&);
    bool _is_correct(const typename Network<FlowType>::EdgeIterator&);
    void _delete_vertex(const size_t&);
    void _change_potential_by_edge(const typename Network<FlowType>::EdgeIterator&,
                                   const bool& push_max = true,
                                   const FlowType& value = FlowType(0));

public:
    MKMMaxFlowFinder(const MKMMaxFlowFinder&) = delete;
    MKMMaxFlowFinder(MKMMaxFlowFinder&&) = delete;
    MKMMaxFlowFinder& operator=(const MKMMaxFlowFinder&) = delete;
    MKMMaxFlowFinder& operator=(MKMMaxFlowFinder&&) = delete;
    explicit MKMMaxFlowFinder(Network<FlowType>&);
};

template<typename FlowType>
MKMMaxFlowFinder<FlowType>::MKMMaxFlowFinder(Network<FlowType>& network):
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
            typename Network<FlowType>::EdgeIterator& edge_it = _direct_iterator[vertex];
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
            typename Network<FlowType>::EdgeIterator& edge_it = _reverse_iterator[vertex];
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
        const typename Network<FlowType>::EdgeIterator& it) {
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
        const typename Network<FlowType>::EdgeIterator& it,
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
    explicit GoldbergMaxFlowFinder(Network<FlowType>&);
    
    //--------------------------------------
    // Your majesty, Find-Max-Flow Function
    //--------------------------------------
    FlowType find_max_flow() override = 0;

protected:
    vector<size_t> _height;
    vector<FlowType> _excess;
    
    bool _is_possible_to_push(const typename Network<FlowType>::EdgeIterator&);
    void _push(const typename Network<FlowType>::EdgeIterator&);
    void _relabel(const size_t&);
    void _push_from_source();
};

template<typename FlowType>
GoldbergMaxFlowFinder<FlowType>::GoldbergMaxFlowFinder(Network<FlowType>& network):
        NetworkMaxFlowFinder<FlowType>(network) {
    _height.resize(network.size(), 0);
    _height[network.get_source()] = network.size();
    _excess.resize(network.size(), FlowType(0));
}


template<typename FlowType>
void GoldbergMaxFlowFinder<FlowType>::_push(
        const typename Network<FlowType>::EdgeIterator& it) {
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
        const typename Network<FlowType>::EdgeIterator& it) {
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
class GoldbergDischargeMaxFlowFinder : protected GoldbergMaxFlowFinder<FlowType> {
public:
    explicit GoldbergDischargeMaxFlowFinder(Network<FlowType>&);
    
    //--------------------------------------
    // Your majesty, Find-Max-Flow Function
    //--------------------------------------
    FlowType find_max_flow() override = 0;

protected:
    bool _discharge(const size_t&);
};

template<typename FlowType>
GoldbergDischargeMaxFlowFinder<FlowType>::GoldbergDischargeMaxFlowFinder(
        Network<FlowType>& network):
        GoldbergMaxFlowFinder<FlowType>(network) {}


template<typename FlowType>
bool GoldbergDischargeMaxFlowFinder<FlowType>::_discharge(const size_t& vertex) {
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
class GoldbergRTFMaxFlowFinder : public GoldbergDischargeMaxFlowFinder<FlowType> {
public:
    explicit GoldbergRTFMaxFlowFinder(Network<FlowType>&);
    
    //--------------------------------------
    // Your majesty, Find-Max-Flow Function
    //--------------------------------------
    FlowType find_max_flow() override;
};

template<typename FlowType>
GoldbergRTFMaxFlowFinder<FlowType>::GoldbergRTFMaxFlowFinder(
        Network<FlowType>& network):
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




template <typename FlowType>
class FordFulkersonMaxFlowFinder : protected NetworkMaxFlowFinder<FlowType> {
private:
    FlowType _answer = FlowType(0);
public:
    explicit FordFulkersonMaxFlowFinder(Network<FlowType>&);
    
    //--------------------------------------
    // Your majesty, Find-Max-Flow Function
    //--------------------------------------
    FlowType find_max_flow() override;
    void bfs(vector<int>&, vector<size_t>&);
};

template<typename FlowType>
FordFulkersonMaxFlowFinder<FlowType>::FordFulkersonMaxFlowFinder(
        Network<FlowType>& network):
            NetworkMaxFlowFinder<FlowType>(network){
    
}


template<typename FlowType>
FlowType FordFulkersonMaxFlowFinder<FlowType>::find_max_flow() {
    vector<int> prev(this->network.size());
    vector<size_t> edges(this->network.size());
    bfs(prev, edges);
    
    while(prev[this->network._target] != -1){
        int u = this->network._target;
        int value = this->network._edges[edges[this->network._target]].get_free_capacity();
        while(u != this->network._source){
            value = std::min(value, this->network._edges[edges[u]].get_free_capacity());
            u = prev[u];
        }
        u = this->network._target;
        while (u != this->network._source) {
            this->network.push(edges[u], value);
            u = prev[u];
        }
        _answer += value;
        bfs(prev, edges);
    }
    return _answer;
}

template<typename FlowType>
void FordFulkersonMaxFlowFinder<FlowType>::bfs(vector<int>& prev, vector<size_t>& edges) {
    prev.assign(prev.size(), -1);
    
    queue<size_t> q;
    
    q.push(this->network._source);
    prev[this->network._source] = -2;
    while (!q.empty()) {
        size_t v = q.front();
        q.pop();
        for (auto it = this->network.begin(v, DIRECT_EDGE);
             it != this->network.end(v, DIRECT_EDGE); ++it) {
            if (prev[it.get_finish()]==-1 && it.get_free_capacity() > FlowType(0)) {
                prev[it.get_finish()] = v;
                edges[it.get_finish()] = *(it.it);
                q.push(it.get_finish());
            }
        }
    }
}


using namespace std;
const int BIGNUMBER = 2500;

#pragma GCC optimize("unroll-loops") // развернуть цикл
#pragma GCC optimize("Ofast") // скомпилировать с о3
#pragma GCC optimize("no-stack-protector") // что-то со стеком

/*
#pragma GCC target("sse,sse2,sse3,ssse3,popcnt,abm,mmx,tune=native") // оптимизации процессора
#pragma GCC optimize("fast-math") // оптимизации сопроцессора
*/
#include <cstring>

void solve() {
    int N, M;
    cin >> N >> M;
    
    int mp[N][M];
    memset(mp, 0, sizeof(mp));
    
    int K, L, x, y;
    cin >> K >> L;
    for (int i = 0; i < K; ++i) {
        cin >> x >> y;
        --x; --y;
        mp[x][y] = 1;
    }
    for (int i = 0; i < L; ++i) {
        cin >> x >> y;
        --x;
        --y;
        mp[x][y] = 2;
    }
    int A_x, A_y, B_x, B_y;
    cin >> A_x >> A_y >> B_x >> B_y;
    --A_x;
    --A_y;
    --B_x;
    --B_y;
    
    Network<int> network(2 * (A_y * N + A_x) + 1, 2 * (B_y * N + B_x), 2 * N * M);
    
    for(int i = 0; i < N; ++i){
        for(int j = 0; j < M; ++j){
            if(i != 0){
                network.add_edge(2 * (j * N + i) + 1, 2 * (j * N + i - 1), BIGNUMBER);
            }
    
            if (j != 0) {
                network.add_edge(2 * (j * N + i) + 1, 2 * ((j - 1) * N + i), BIGNUMBER);
            }
    
            if (i != N - 1) {
                network.add_edge(2 * (j * N + i) + 1, 2 * (j * N + i + 1), BIGNUMBER);
            }
    
            if (j != M - 1) {
                network.add_edge(2 * (j * N + i) + 1, 2 * ((j + 1) * N + i), BIGNUMBER);
            }
            if(mp[i][j] == 0){
                network.add_edge(2 * (j * N + i), 2 * (j * N + i) + 1, BIGNUMBER);
            }
            if (mp[i][j] == 2) {
                network.add_edge(2 * (j * N + i), 2 * (j * N + i) + 1, 1);
            }
        }
    }
    
    FordFulkersonMaxFlowFinder<int> solver(network);
    
    
    int answer = solver.find_max_flow();
    
    if(answer >= BIGNUMBER){
        cout << -1;
        return;
    }
    
    cout << answer << endl;
    vector<bool> ans;
    network.find_reachable(ans);
    
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < M; ++j) {
            if((i != A_x || j != A_y) && (i != B_x || j != B_y) &&
                ans[2*(j * N + i)] && !ans[2 * (j * N + i) + 1] && mp[i][j] == 2){
                cout << i + 1 << ' ' << j + 1 << endl;
            }
        }
    }
    
    /*
#ifdef LOCAL
    for (auto i : network._edges) {
        cout << i.get_start() + 1 << "->" << i.get_finish() + 1 << ' ' << i._flow << '/' << i._capacity << endl;
    }
#endif*/
}

void test(){
    Network<int> network(0, 3, 4);
    network.add_edge(0, 1, 5);
    network.add_edge(0, 2, 4);
    network.add_edge(1, 2, 1);
    network.add_edge(1, 3, 3);
    network.add_edge(2, 3, 7);
    
    FordFulkersonMaxFlowFinder<int> solver(network);
    
    int answer = solver.find_max_flow();
    
    cout << answer << endl;
}

int main() {
#ifdef LOCAL
    freopen("input.txt", "rt", stdin);
    //freopen("output.txt", "wt", stdout);
#endif
    cin.tie(0);
    ios_base::sync_with_stdio(0);
    solve();
    //test();
}