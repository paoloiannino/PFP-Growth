
***think about confidence/correlation***
# Parallel implementation of FP Growth Algorithm

This project is part of the course *High Performance Computing (2017)* at ENSIMAG (Grenoble).

## Getting Started

This section will contain the instruction to evaluate the project.


## FP Growth Pseudo Code

```
Input: DB (Containing transaction), N
Output: Top N frequent itemset

- build the occurrency vector/list and order the element with respect to the support count
- order the items in the transactions (in the DB) with respect to the order obtained by the previous step
- create the tree having the root equal to null element and iterating over the transaction on the ordered DB
- compute the frequent itemset starting from the leaf and counting the supports

```

## Authors

* **Iannino Paolo**
* **La Quatra Moreno**

## License

See the [LICENSE](LICENSE) file for details
