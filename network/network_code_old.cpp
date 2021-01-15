#include <iostream>
#include <vector>
#include <queue>
#include <list>


using std::vector;
using std::queue;
using std::list;


template<typename FlowType>
class Network{
public:
    // Network own methods

    Network() = delete;
    explicit Network(const size_t&, const size_t&, const size_t&);

    [[nodiscard]] const size_t& get_source() const;
    [[nodiscard]] const size_t& get_target() const;
    size_t size();

    //--------------------------------------
    // Your majesty, Find-Max-Flow Function
    //--------------------------------------
    template<template<typename> class SearchEngine>
    FlowType find_max_flow();


    // Methods connected with Edge

    class Edge;
    void add_edge(const size_t&, const size_t&, const FlowType&);
    const Edge& get_edge(const size_t&) const;
    void push(const size_t&, const FlowType&);


    // Methods connected with EdgeIterator

    enum EDGE_ITERATOR_TYPE{
        DIRECT_GRAPH_ITERATOR = false,
        REVERSE_GRAPH_ITERATOR = true
    };

    class EdgeIterator;
    EdgeIterator begin(const size_t&, EDGE_ITERATOR_TYPE);
    EdgeIterator end(const size_t&, EDGE_ITERATOR_TYPE);

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
class Network<FlowType>::Edge{
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
    // Network<FlowType>& _network;
};

template<typename FlowType>
Network<FlowType>::Edge::Edge(const size_t& start, const size_t& finish, const FlowType& capacity): _start(start),
                                                                                                    _finish(finish),
                                                                                                    _capacity(capacity){}

template<typename FlowType>
void Network<FlowType>::Edge::push(const FlowType& value) {
    // std::cout << "[Flow updated from " << _flow << " to " << _flow + value << ", capacity = " << _capacity << "]\n";
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
class Network<FlowType>::EdgeIterator: public std::iterator<std::forward_iterator_tag,
                                                            typename Network<FlowType>::Edge>{
public:
    EdgeIterator() = delete;
    explicit EdgeIterator(const size_t&, Network<FlowType>&, const bool&);
    explicit EdgeIterator(const vector<size_t>::iterator&, Network<FlowType>&);

    bool operator!=(EdgeIterator const& other) const;
    bool operator==(EdgeIterator const& other) const;

    EdgeIterator& operator++();

    FlowType get_free_capacity() const;
    [[nodiscard]] const size_t& get_start() const;
    [[nodiscard]] const size_t& get_finish() const;

    void push(const FlowType&) const;
    void push_max() const;
private:
    vector<size_t>::iterator it;
    Network& network;
};

template<typename FlowType>
Network<FlowType>::EdgeIterator::EdgeIterator(const size_t& vertex_index, Network<FlowType>& network, const bool& reverse):
    it((reverse ? network._rev_graph : network._graph)[vertex_index].begin()),
    network(network){}

template<typename FlowType>
Network<FlowType>::EdgeIterator::EdgeIterator(const vector<size_t>::iterator& it, Network& network):
    it(it), network(network) {}

template<typename FlowType>
bool Network<FlowType>::EdgeIterator::operator==(const EdgeIterator& other) const {
    return it == other.it;
}

template <typename FlowType>
bool Network<FlowType>::EdgeIterator::operator!=(const EdgeIterator& other) const {
    return it != other.it;
}

template<typename FlowType>
typename Network<FlowType>::EdgeIterator& Network<FlowType>::EdgeIterator::operator++() {
    ++it;
    return *this;
}

template<typename FlowType>
FlowType Network<FlowType>::EdgeIterator::get_free_capacity() const {
    return network.get_edge(*it).get_free_capacity();
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
Network<FlowType>::Network(const size_t& source, const size_t& target, const size_t& vertex_number):
    _source(source),
    _target(target) {
    int a;
    // ToDo: add exceptions for situations when source|target >= vertex_number
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

//--------------------------------------
// Your majesty, Find-Max-Flow Function
//--------------------------------------
template<typename FlowType>
template<template<typename> class SearchEngine>
FlowType Network<FlowType>::find_max_flow() {
    SearchEngine<FlowType> search_engine(*this);
    return search_engine.solve();
}


template<typename FlowType>
void Network<FlowType>::add_edge(const size_t& start, const size_t& finish, const FlowType& capacity) {
    _edges.emplace_back(start, finish, capacity);
    _edges.emplace_back(finish, start, 0);

    _graph[start].push_back(_edges.size() - 2);
    _graph[finish].push_back(_edges.size() - 1);

    _rev_graph[finish].push_back(_edges.size() - 2);
    _rev_graph[start].push_back(_edges.size() - 1);
}

template<typename FlowType>
const typename Network<FlowType>::Edge& Network<FlowType>::get_edge(const size_t& edge_index) const {
    return _edges[edge_index];
}

template<typename FlowType>
void Network<FlowType>::push(const size_t& edge_index, const FlowType& value) {
    _edges[edge_index].push(value);
    _edges[edge_index ^ 1u].push(-value);
}



template<typename FlowType>
typename Network<FlowType>::EdgeIterator Network<FlowType>::begin(const size_t& vertex_index,
                                                                  Network<FlowType>::EDGE_ITERATOR_TYPE type) {
    return EdgeIterator(vertex_index, *this, type);
}

template<typename FlowType>
typename Network<FlowType>::EdgeIterator Network<FlowType>::end(const size_t& vertex_index,
                                                                  Network<FlowType>::EDGE_ITERATOR_TYPE type) {
    return EdgeIterator((type ? _rev_graph : _graph)[vertex_index].end(), *this);
}














// -----------------------------------------
//           Network Search Engines
// -----------------------------------------

template<typename FlowType>
class MaxFlowSearchEngine{
protected:
    Network<FlowType>& network;
    size_t source;
    size_t target;
public:
    explicit MaxFlowSearchEngine(Network<FlowType>&);
    virtual FlowType solve() = 0;
};

template<typename FlowType>
MaxFlowSearchEngine<FlowType>::MaxFlowSearchEngine(Network<FlowType>& network): network(network),
                                                                                source(network.get_source()),
                                                                                target(network.get_target()){}


// -----------------------------------------------------
//   Malhotra-Kumar-Maheshwari Network Search Engine
// -----------------------------------------------------
/* *
 * Find maximum flow in arbitary network for O(V^3)
 * Using concept of blocking flows, it finds blocking flow
 * in O(V^2) and make O(V) iterations of searching blocking
 * flow (O(V^3) as a whole)
 */

template<typename FlowType>
class MKMMaxFlowSE : protected MaxFlowSearchEngine<FlowType>{
    typedef MaxFlowSearchEngine<FlowType> SE; // shortcut for SearchEngine to get
                                                // shorter expression to access network field
private:
    vector<FlowType> _potential_in;
    vector<FlowType> _potential_out;
    vector<bool> _deleted;
    vector<size_t> _level;
    vector<typename Network<FlowType>::EdgeIterator> _direct_iterator, _reverse_iterator;
    FlowType answer = FlowType(0);

    void _calculate_level();
    void _calculate_potential();
    void _get_potential(const size_t&);
    size_t _find_min_potential_vertex();
    void _find_blocking_flow();
    void _prepare_for_new_iter();
    void _push_to_target(const size_t&);
    void _push_to_source(const size_t&);
    void _is_correct(const typename Network<FlowType>::EdgeIterator&);
    void _delete_vertex(const size_t&);
    void _change_potential_by_edge(const typename Network<FlowType>::EdgeIterator&,
                                   const bool& push_max = true,
                                   const FlowType& value = FlowType(0));

public:
    MKMMaxFlowSE(const MKMMaxFlowSE&) = delete;
    MKMMaxFlowSE(MKMMaxFlowSE&&) = delete;
    MKMMaxFlowSE& operator= (const MKMMaxFlowSE&) = delete;
    MKMMaxFlowSE& operator= (MKMMaxFlowSE&&) = delete;

    explicit MKMMaxFlowSE(Network<FlowType>&);

    FlowType solve() override;
};

template<typename FlowType>
MKMMaxFlowSE<FlowType>::MKMMaxFlowSE(Network<FlowType>& network): MaxFlowSearchEngine<FlowType>(network){
    _potential_in.resize(network.size(), FlowType(0));
    _potential_out.resize(network.size(), FlowType(0));
    _deleted.resize(network.size(), false);
    _level.resize(network.size(), 0);
    for(size_t i = 0; i < SE::network.size(); ++i){
        _direct_iterator.push_back(SE::network.begin(i, Network<FlowType>::DIRECT_GRAPH_ITERATOR));
        _reverse_iterator.push_back(SE::network.begin(i, Network<FlowType>::REVERSE_GRAPH_ITERATOR));
    }
}

template<typename FlowType>
void MKMMaxFlowSE<FlowType>::_calculate_level() {
    vector<bool> visited(SE::network.size());
    queue<size_t> to_visit;
    to_visit.push(SE::source);
    _level[SE::source] = 0;

    while(!to_visit.empty()) {
        size_t vertex = to_visit.front();
        to_visit.pop();
        if(!visited[vertex]) {
            visited[vertex] = true;
            for(auto vertex_edge_it = SE::network.begin(vertex, Network<FlowType>::DIRECT_GRAPH_ITERATOR);
                vertex_edge_it != SE::network.end(vertex);
                ++vertex_edge_it){
                size_t finish = vertex_edge_it.get_finish();
                if(!visited[finish]) {
                    _level[finish] = _level[vertex] + 1;
                    to_visit.push(finish);
                }
            }
        }
    }
}

template<typename FlowType>
void MKMMaxFlowSE<FlowType>::_calculate_potential() {
    for(size_t vertex = 0; vertex < SE::network.size(); ++vertex){
        for(auto vertex_edge_it = SE::network.begin(vertex, Network<FlowType>::DIRECT_GRAPH_ITERATOR);
            vertex_edge_it != SE::network.end(vertex);
            ++vertex_edge_it){
            if(_is_correct(vertex_edge_it)){
                FlowType capacity = vertex_edge_it.get_free_capacity();
                _potential_out[vertex] += capacity;
                _potential_in[vertex_edge_it.get_finish()] += capacity;
            }
        }
    }
}


template<typename FlowType>
void MKMMaxFlowSE<FlowType>::_get_potential(const size_t& vertex) {
    if(vertex == SE::source) {
        return _potential_out[vertex];
    }
    if(vertex == SE::source) {
        return _potential_in[vertex];
    }
    return std::min(_potential_in[vertex], _potential_out[vertex]);
}

template<typename FlowType>
size_t MKMMaxFlowSE<FlowType>::_find_min_potential_vertex() {
    size_t min_potential_vertex = SE::source;
    for(size_t i = 0; i < SE::network.size(); ++i){
        if(!_deleted[i] && _get_potential(i) < _get_potential(min_potential_vertex)){
            min_potential_vertex = i;
        }
    }
    return min_potential_vertex;
}

template<typename FlowType>
void MKMMaxFlowSE<FlowType>::_find_blocking_flow() {
    // maybe !_deleted[SE::source] && !_deleted[SE::network.get_target()])
    while(!_deleted[SE::source]) {
        size_t min_potential_vertex = _find_min_potential_vertex();
        if(_get_potential(min_potential_vertex) == FlowType(0)){
            _delete_vertex(min_potential_vertex);
        } else {
            answer += _get_potential(min_potential_vertex);
            _push_to_source(min_potential_vertex);
            _push_to_target(min_potential_vertex);
        }
    }
}

template<typename FlowType>
void MKMMaxFlowSE<FlowType>::_prepare_for_new_iter() {
    _deleted.assign(_deleted.size(), false);
    _potential_in.assign(_potential_in.size(), 0);
    _potential_out.assign(_potential_out.size(), 0);
    for(size_t i = 0; i < SE::network.size(); ++i){
        _direct_iterator[i] = SE::network.begin(i, Network<FlowType>::DIRECT_GRAPH_ITERATOR);
        _reverse_iterator[i] = SE::network.begin(i, Network<FlowType>::REVERSE_GRAPH_ITERATOR);
    }
    _calculate_level();
    _calculate_potential();
}

template<typename FlowType>
void MKMMaxFlowSE<FlowType>::_push_to_target(const size_t& min_potential_vertex) {
    vector<FlowType> to_be_pushed(SE::network.size());
    queue<size_t> blocking_flow_vertexes;
    blocking_flow_vertexes.push(min_potential_vertex);
    to_be_pushed[min_potential_vertex] = _get_potential(min_potential_vertex);


    while(!blocking_flow_vertexes.empty()){
        size_t vertex = blocking_flow_vertexes.front();
        blocking_flow_vertexes.pop();
        while(to_be_pushed[vertex] != FlowType(0)){
            typename Network<FlowType>::EdgeIterator& edge_it = _direct_iterator[vertex];
            if(_is_correct(edge_it)){
                if(to_be_pushed[edge_it.get_finish()] == FlowType(0)){
                    blocking_flow_vertexes.push(edge_it.get_finish());
                }
                if(edge_it.get_free_capacity() > to_be_pushed[vertex]){
                    to_be_pushed[edge_it.get_finish()] += to_be_pushed[vertex];
                    edge_it.push(to_be_pushed[vertex]);
                    to_be_pushed[vertex] = 0;
                    break;
                }
                to_be_pushed[edge_it.get_finish()] += edge_it.get_free_capacity();
                to_be_pushed[vertex] -= edge_it.get_free_capacity();
                edge_it.push_max();
            }
            ++_direct_iterator[vertex];
        }
    }
}

template<typename FlowType>
void MKMMaxFlowSE<FlowType>::_push_to_source(const size_t& min_potential_vertex) {
    vector<FlowType> to_be_pushed(SE::network.size());
    queue<size_t> blocking_flow_vertexes;
    blocking_flow_vertexes.push(min_potential_vertex);
    to_be_pushed[min_potential_vertex] = _get_potential(min_potential_vertex);


    while(!blocking_flow_vertexes.empty()){
        size_t vertex = blocking_flow_vertexes.front();
        blocking_flow_vertexes.pop();
        while(to_be_pushed[vertex] != FlowType(0)){
            typename Network<FlowType>::EdgeIterator& edge_it = _direct_iterator[vertex];
            if(_is_correct(edge_it)){
                if(to_be_pushed[edge_it.get_finish()] == FlowType(0)){
                    blocking_flow_vertexes.push(edge_it.get_finish());
                }
                if(edge_it.get_free_capacity() > to_be_pushed[vertex]){
                    _change_potential_by_edge(edge_it, false, to_be_pushed[vertex]);
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
void MKMMaxFlowSE<FlowType>::_is_correct(const typename Network<FlowType>::EdgeIterator& it) {
    return _level[it.get_start()] + 1 == _level[it.get_finish()] &&
           !_deleted[it.get_start()] &&
           !_deleted[it.get_finish()] &&
           it.get_free_capacity() > FlowType(0);
}


template<typename FlowType>
void MKMMaxFlowSE<FlowType>::_delete_vertex(const size_t& vertex) { // delete if potential is 0
    if(_get_potential(vertex) != FlowType(0)){
        return;
    }
    auto iterator_type = _potential_in[vertex] == FlowType(0) ?
                         Network<FlowType>::DIRECT_GRAPH_ITERATOR :
                         Network<FlowType>::REVERSE_GRAPH_ITERATOR;


    for (auto vertex_edge_it = SE::network.begin(vertex, iterator_type);
         vertex_edge_it != SE::network.end(vertex);
         ++vertex_edge_it){
        if(!_is_correct(vertex_edge_it)) {
            _change_potential_by_edge(vertex_edge_it);
        }
    }
    _deleted[vertex] = true;
}


template<typename FlowType>
FlowType MKMMaxFlowSE<FlowType>::solve() {
    _prepare_for_new_iter();
    while(_level[SE::network.start]){
        _find_blocking_flow();
        _prepare_for_new_iter();
    }
    return answer;
}

template<typename FlowType>
void MKMMaxFlowSE<FlowType>::_change_potential_by_edge(const typename Network<FlowType>::EdgeIterator& it,
                                                       const bool& push_max,
                                                       const FlowType& value) {
    if(push_max){
        _potential_in[it.get_finish()] -= it.get_free_capacity();
        _potential_out[it.get_start()] -= it.get_free_capacity();
    } else {
        _potential_in[it.get_finish()] -= value;
        _potential_out[it.get_start()] -= value;
    }
}









// ---------------------------------------------------------
//      Goldberg [push-relable] Network Search Engine
// ---------------------------------------------------------

/* *
 * Alternative Goldberg technology for transforming pre-flow
 * into max flow. Standard algorithm works in O(EV^2). Its
 * improved version [relabel-to-front, discharge] technology
 * works in O(V^3)
 */

template<typename FlowType>
class GoldbergMaxFlowSE : protected MaxFlowSearchEngine<FlowType>{
    typedef MaxFlowSearchEngine<FlowType> SE; // shortcut for SearchEngine to get
                                                // shorter expression to access network field
private:
    vector<size_t> height;
    vector<FlowType> excess;

    bool _is_possible_to_push(const typename Network<FlowType>::EdgeIterator&);
    void _push(const typename Network<FlowType>::EdgeIterator&);
    void _relabel(const size_t&);
    bool _discharge(const size_t&);
    void _push_from_source();

public:
    explicit GoldbergMaxFlowSE(Network<FlowType>&);
    FlowType solve() override;
};

template<typename FlowType>
GoldbergMaxFlowSE<FlowType>::GoldbergMaxFlowSE(Network<FlowType>& network): MaxFlowSearchEngine<FlowType>(network) {
    height.resize(network.size(), 0);
    height[network.get_source()] = network.size();
    excess.resize(network.size(), FlowType(0));
}

template<typename FlowType>
void GoldbergMaxFlowSE<FlowType>::_push(const typename Network<FlowType>::EdgeIterator& it) {
    FlowType value = std::min(it.get_free_capacity(), excess[it.get_start()]);
    it.push(value);
    excess[it.get_finish()] += value;
    excess[it.get_start()] -= value;
}

template<typename FlowType>
void GoldbergMaxFlowSE<FlowType>::_relabel(const size_t& vertex) {
    size_t new_height = 2 * SE::network.size() + 1;
    for (auto it = SE::network.begin(vertex, Network<FlowType>::DIRECT_GRAPH_ITERATOR);
         it != SE::network.end(vertex, Network<FlowType>::DIRECT_GRAPH_ITERATOR);
         ++it) {
        if(it.get_free_capacity() > FlowType(0)){
            new_height = std::min(new_height, height[it.get_finish()] + 1);
        }
    }
    height[vertex] = new_height;
}

template<typename FlowType>
bool GoldbergMaxFlowSE<FlowType>::_discharge(const size_t& vertex) {
    bool discharged = false;
    for (auto it = SE::network.begin(vertex, Network<FlowType>::DIRECT_GRAPH_ITERATOR);
         it != SE::network.end(vertex, Network<FlowType>::DIRECT_GRAPH_ITERATOR);
         ++it) {
        if(_is_possible_to_push(it)){
            _push(it);
            if(excess[vertex] <= FlowType(0)){
                discharged = true;
                break;
            }
        }
    }
    return discharged;
}

template<typename FlowType>
bool GoldbergMaxFlowSE<FlowType>::_is_possible_to_push(const typename Network<FlowType>::EdgeIterator& it) {
    return (height[it.get_finish()] + 1 == height[it.get_start()]) && (it.get_free_capacity() > FlowType(0));
}

template<typename FlowType>
void GoldbergMaxFlowSE<FlowType>::_push_from_source() {
    size_t source = SE::source;
    FlowType source_excess = FlowType(0);
    for (auto it = SE::network.begin(source, Network<FlowType>::DIRECT_GRAPH_ITERATOR);
         it != SE::network.end(source, Network<FlowType>::DIRECT_GRAPH_ITERATOR);
         ++it){
        source_excess += it.get_free_capacity();
        excess[source] += it.get_free_capacity();
        _push(it);
    }
    excess[source] -= source_excess;
}



template<typename FlowType>
FlowType GoldbergMaxFlowSE<FlowType>::solve() {
    _push_from_source();

    list<size_t> vertexes;
    for(size_t i = 0; i < SE::network.size(); ++i){
        if(i != SE::source && i != SE::network.get_target()) {
            vertexes.push_back(i);
        }
    }

    auto it = vertexes.begin();
    while(it != vertexes.end()){
        //std::cout << *it << std::endl;
        if((excess[*it] == FlowType(0)) || _discharge(*it)){
            ++it;
        } else {
            _relabel(*it);
            size_t moved = *it;
            vertexes.erase(it);
            vertexes.push_front(moved);
            it = vertexes.begin();
        }
    }

    return -excess[SE::source];
}


using std::cout;
using std::cin;
using std::endl;

void solve(){
    size_t n, m;
    cin >> n;
    Network<int> network(0, n + 1, n + 2);
    int profit, sum_profit = 0;
    for(size_t i = 0; i < n; ++i){
        cin >> profit;
        if(profit >= 0){
            network.add_edge(0, i + 1, profit);
            sum_profit += profit;
        }
        else {
            network.add_edge(i + 1, n + 1, -profit);
        }
    }
    size_t dependence;
    for(size_t i = 0; i < n; ++i){
        cin >> m;
        for(size_t j = 0; j < m; ++j){
            cin >> dependence;
            network.add_edge(i + 1, dependence, 1000000);
        }
    }

    //auto pre_ans = network.find_max_flow<MKMMaxFlowSE>();
    auto pre_ans = network.find_max_flow<GoldbergMaxFlowSE>();
    /*cout << pre_ans << endl;
    cout << "\nAnswer is " << sum_profit - pre_ans << endl;*/

    cout << sum_profit - pre_ans << endl;
}

int main() {
#ifdef LOCAL
    freopen("input.txt", "rt", stdin);
    //freopen("output.txt", "wt", stdout);
#endif
// something
    solve();
}
