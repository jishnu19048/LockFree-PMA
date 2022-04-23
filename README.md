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

- [Jishnu Raj Parashar](https://github.com/jishnu19048) <br>
- [Nandika Jain](https://github.com/nandikajain) <br>
- [Saatvik Bhatnagar](https://github.com/Saatvik07) 

### References
*[1]* M.A. Bender, E.D. Demaine, and M. Farach-Colton. Cache-oblivious b-trees. In Proceedings 41st Annual Symposium on Foundations of Computer Science, pages 399–409, 2000.<br>
*[2]* Michael A. Bender, Jeremy T. Fineman, Seth Gilbert, and Bradley C. Kuszmaul. Concurrent cache-oblivious b-trees. In Proceedings of the Seventeenth Annual ACM Symposium on Parallelism in Algorithms and Architectures, SPAA ’05, page 228–237, New York, NY, USA, 2005. Association for Computing Machinery.
