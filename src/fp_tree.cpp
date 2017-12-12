#include "fp_tree.h"

fp_tree::fp_tree()
{
    this->root_node = make_shared<fp_node>();
    this->root_node->parent = NULL;
    this->root_node->support = -1;
    this->root_node->item = string();
    this->pointer_table = make_shared<map<string, vector<shared_ptr<fp_node>>>>();
}

fp_tree::~fp_tree()
{
    clean_tree(root_node);
}

void fp_tree::insert_transaction(vector<string>::const_iterator next_item, vector<string>::const_iterator transaction_end)
{
    this->insert(this->root_node, next_item, transaction_end);
}

void fp_tree::clean_tree(shared_ptr<fp_node> node){
    if(node->childs.size() == 0){
        node->childs.clear();
        return;
    }

    for(shared_ptr<fp_node> child : node->childs){
        clean_tree(child);
        child.reset();
    }

    node->childs.clear();
    return;

}

shared_ptr<vector<vector<string>>> fp_tree::get_transaction(string item)
{
    vector<shared_ptr<fp_node>> pointer_vector = pointer_table->operator[](item);
    shared_ptr<vector<vector<string>>> transactions = make_shared<vector<vector<string>>>();
    pair<shared_ptr<vector<string>>, int> transaction_support;

    for(vector<shared_ptr<fp_node>>::const_iterator it = pointer_vector.begin(); it != pointer_vector.end(); it++){
        transaction_support = this->get(*it);
        for(int i = 0; i < transaction_support.second; i++)
            transactions->push_back(*transaction_support.first);
    }


    return transactions;
}


void fp_tree::insert(shared_ptr<fp_node> current_node, vector<string>::const_iterator next_item, vector<string>::const_iterator transaction_end)
{
    shared_ptr<fp_node> next_node = NULL;

    if(next_item == transaction_end)
        return;

    for(shared_ptr<fp_node> child : current_node->childs){
        if(next_item->compare(child->item) == 0){
            next_node = child;
            break;
        }
    }

    if(next_node == NULL){
        next_node = make_shared<fp_node>();
        next_node->item = *next_item;
        next_node->support = 1;
        next_node->parent = current_node;
        current_node->childs.push_back(next_node);
        pointer_table->operator[](*next_item).push_back(next_node);
    } else {
        next_node->support++;
    }

    fp_tree::insert(next_node, ++next_item, transaction_end);
}

pair<shared_ptr<vector<string>>, int> fp_tree::get(shared_ptr<fp_node> start_node)
{
    shared_ptr<vector<string>> transaction = make_shared<vector<string>>();

    for(shared_ptr<fp_node> node_pointer = start_node->parent; node_pointer->parent != NULL; node_pointer = node_pointer->parent)
           transaction->push_back(node_pointer->item);

    return pair<shared_ptr<vector<string>>, int>(transaction, start_node->support);
}


