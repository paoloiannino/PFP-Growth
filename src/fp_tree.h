#ifndef FP_TREE_H
#define FP_TREE_H

#include <istream>
#include <vector>
#include <memory>
#include <map>

using namespace std;

struct fp_node {
    string item;
    int support;
    vector<shared_ptr<fp_node>> childs;
    shared_ptr<fp_node> parent;
};

class fp_tree
{
    public:
	fp_tree();
    ~fp_tree();
    void insert_transaction(vector<string>::const_iterator next_item, vector<string>::const_iterator transaction_end);
    shared_ptr<vector<vector<string>>> get_transaction(string item);

    private:
    shared_ptr<fp_node> root_node;
    shared_ptr<map<string, vector<shared_ptr<fp_node>>>> pointer_table;

    void insert(shared_ptr<fp_node> current_node, vector<string>::const_iterator next_item, vector<string>::const_iterator transaction_end);
    pair<shared_ptr<vector<string>>, int> get(shared_ptr<fp_node> start_node);
    void clean_tree(shared_ptr<fp_node> node);
};

#endif // FP_TREE_H
