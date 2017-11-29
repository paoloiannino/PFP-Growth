**PseudoCode FP Growth**
***think about confidence/correlation***


INPUT: transactions, pattern, threshold
*transactions:  is a Vector of Vectors of strings (each Row is a Transaction, each string is an item)*
*pattern:       is a Vector of strings that WILL contain the conditional pattern*
*threshold:     contains the minimum support required for a frequent itemset*

OUTPUT: frequent_itemset
*frequent_itemset:   is a vector of c++ pairs containing < transactions (vector of strings), support (int) >*
