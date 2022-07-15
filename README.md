# SortingAlgorithmVisualizer
Includes a module with some sorting algorithms and the main script allows to select a one of the algorithms and visualize the process of sorting.

Author: David García Allo
Date: 11/07/2022

Description: 
Showing with bars the process followed by the different sorting algorithms.
Choose an algorithm contained in SortingAlgorithms.py to sort the given array and visualize each step.
Press 'space' to init the sorting process
Press 'c' while not executing the algorithm to create a new random array.
Module 'pygame' needed, if needed use command 'pip install pygame'

Contents:
- README.md: General information.
- SortingAlgorithms.py: Module with the different algos that we can use:
    {Bubble, Selection, Insertion, Merge, Quicksort, Counting, Radix, Bucket, Heap, Shell}
- ascending.ico: Icon of the pygame window
- main.py: Main script that contains all the functions, classes and modules needed to visualize the sorting process

Ideas to improve:
- Maybe add some buttons that allows to execute algorithm, reset array and kill the process
- Improve the visualixation of the sorting methods that uses different data structures that doesn't allow us a proper visualization.
- Add some check-box or similar to choose one of the available algos
