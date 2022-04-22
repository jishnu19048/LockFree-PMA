## Lock Free PMA
We extend the existing sequential PMA to a parallel lock-free setting by using compare and swap (CAS). The data structure supports lock-free searching, insertions, and deletions. The data structure is linearizable, i.e., we can express a set of operations in a parallel setting in a serialized order. The data structure preserves the set of operations as the original sequential PMA, and gives a better performance with parallelism.

### Running Instructions
> 1. Clone this repository using the command: `git clone https://github.com/jishnu19048/LockFree-PMA.git`
> 2. cd into the directory by ```cd LockFree-PMA```.
> 3. The concurrent code is in the branch ````concurrent````, checkout the branch by ```git checkout concurrent```.
> 4. Run ```gcc pma_v1.c -o run``` and then ```./run``` to run the code.

### Dependencies
1. A GNU Compiler Collection to compile the C program.
2. Git configured for easy access and branch switching.

### *Contact us:*

- [Jishnu Raj Parashar Jain](https://github.com/jishnu19048) <br>
- [Nandika Jain](https://github.com/nandikajain) <br>
- [Saatvik Bhatnagar](https://github.com/Saatvik07) 
