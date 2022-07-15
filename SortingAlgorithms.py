"""
Different Sorting Algorithms

Bubble, Selection, Insertion, Merge, Quicksort, Counting, Radix, Bucket, Heap,
Shell (Different orders of complexity and stability)

Reference: https://www.programiz.com/dsa/sorting-algorithm

David García Allo - 25/06/2022
"""

def BubbleSort(array, draw_array = None):
    """              Bubble sort: O(n·n); Stable
    
    Bubble sort is a sorting algorithm that compares two adjacent elements
    and swaps them until they are not in the intended order.
    
    SCHEME:
    bubbleSort(array)
     for i <- 1 to indexOfLastUnsortedElement-1
      if leftElement > rightElement
       swap leftElement and rightElement
    end bubbleSort
    
    APPLICATIONS:
    Bubble sort is used if
      complexity does not matter
      short and simple code is preferred"""
      
    size = len(array)
    #Loop through all elements
    for i in range(size): 
        swap = False
        #Loop for comparison
        for j in range(size - 1 - i): 
            #Swaping elements, bigger element to the right
            if array[j] > array[j+1]: 
                temp = array[j]
                array[j] = array[j+1]
                array[j+1] = temp
                swap = True
            if draw_array:
                draw_array()
        #If the last iteration didn't swap any element, it is sorted
        if not swap: 
            return array
        
    return array



def SelectionSort(array, draw_array = None):
    """              Selection sort: O(n·n); Non-Stable
    
    Selection sort is a sorting algorithm that selects the smallest element 
    from an unsorted list in each iteration and places that element at the 
    beginning of the unsorted list.
    
    SCHEME:
    selectionSort(array)
     repeat (size - 1) times
     set the first unsorted element as the minimum
     for each of the unsorted elements
      if element < currentMinimum
       set element as new minimum
     swap minimum with first unsorted position
    end selectionSort
    
    APPLICATIONS
    The selection sort is used when
      a small list is to be sorted
      cost of swapping does not matter
      checking of all the elements is compulsory
      cost of writing to a memory matters like in flash memory (number of 
        writes/swaps is O(n) as compared to O(n2) of bubble sort)"""
      
      
    size = len(array)
    #Loop through all elements
    for i in range(size):
        minimum = i #Set the first unsorted element as the minimum
        #Loop through unsorted elements
        for j in range(i+1, size):
            #Setting new minimum
            if array[j]<array[minimum]:
                minimum = j
            if draw_array:
                draw_array()
        #Swapping minimum to the first unsorted element
        temp = array[i]
        array[i] = array[minimum]
        array[minimum] = temp
        if draw_array:
                draw_array()
    return array


def InsertionSort(array, draw_array=None):
    """              Insertion sort: O(n·n); Stable
    
    Insertion sort is a sorting algorithm that places an unsorted element at
    its suitable place in each iteration. Insertion sort works similarly as we 
    sort cards in our hand in a card game.
    
    SCHEME:
    
    insertionSort(array)
     mark first element as sorted
     for each unsorted element X
      'extract' the element X
      for j <- lastSortedIndex down to 0
       if current element j > X
        move sorted element to the right by 1
      break loop and insert X here
    end insertionSort
    
    APPLICATIONS:
    The insertion sort is used when:
        the array is has a small number of elements
        there are only a few elements left to be sorted"""
        
    size = len(array)
    #Loop through all elements, except the first
    for i in range(1, size):
        unsortElement = array[i] #Extracting the element
        j = i - 1
        
        #While loop better than reversed for loop because it allows reach j=-1
        #which it is important to place the lowest value at array[0]
        while j>=0 and unsortElement < array[j]:
            array[j+1] = array[j]
            j = j - 1
            if draw_array:
                draw_array()
                
        array[j+1] = unsortElement #Insrtion in the correct place
        if draw_array:
                draw_array()
        
    return array



def MergeSort(array, draw_array = None):
    """              Merge sort: O(n·log n); Stable
    
    Merge Sort is one of the most popular sorting algorithms that is based 
    on the principle of Divide and Conquer Algorithm. The merge sort algorithm 
    recursively divides the array into halves until we reach the base case of 
    array with 1 element. After that, the merge function picks up the sorted 
    sub-arrays and merges them to gradually sort the entire array.
    
    SCHEME:
    MergeSort(A, p, r):
      if p > r 
        return
      q = (p+r)/2
      mergeSort(A, p, q)
      mergeSort(A, q+1, r)
      merge(A, p, q, r)
      
    Merge(A, p, q, r):
        Have we reached the end of any of the arrays?
    No:
        Compare current elements of both arrays 
        Copy smaller element into sorted array
        Move pointer of element containing smaller element
    Yes:
        Copy all remaining elements of non-empty array
        
    APPLICATIONS:
    Inversion count problem
    External sorting
    E-commerce applications"""
    
    #If len = 1 mergeSort will do nothing, it's the stop of the division
    if len(array) > 1:
        #r is the element which we use to divide arrays
        r = len(array)//2
        A1 = array[:r]
        A2 = array[r:]
        
        #Sort each subarray
        A1 = MergeSort(A1)
        A2 = MergeSort(A2)
        #We will only pass this lines when len = 1
        
        #MERGE
        #Index for the merge step
        i = j = k = 0
        #Compare elements of both arrays
        while i < len(A1) and j < len(A2):
            if A1[i] < A2[j]:
                array[k] = A1[i]
                i += 1
                k += 1
                if draw_array:
                    draw_array()
            else:
                array[k] = A2[j]
                j += 1
                k += 1
                if draw_array:
                    draw_array()
            
        
        #When we get to the end of one of the arrays, just copy the other one
        #which it's already sorted
        while i < len(A1):
            array[k] = A1[i]
            i += 1
            k += 1
            if draw_array:
                draw_array()
            
        while j < len(A2):
            array[k] = A2[j]
            j += 1
            k += 1
            if draw_array:
                draw_array()
        
        return array #Merge of the sorted subarrays

    return array #If len=1 -> only 1 element -> Already sorted


def QuickSort(array, lowestIndex = 0, highestIndex = False, draw_array = None):
    """              Quicksort: O(n·log n); Non-Stable
    
    Quicksort is a sorting algorithm based on the divide and conquer approach 
    where: 
        An array is divided into subarrays by selecting a pivot element 
        (element selected from the array).
        
        While dividing the array, the pivot element should be positioned in 
        such a way that elements less than pivot are kept on the left side and 
        elements greater than pivot are on the right side of the pivot.
        
        The left and right subarrays are also divided using the same approach. 
        This process continues until each subarray contains a single element.
        
        At this point, elements are already sorted. Finally, elements are
        combined to form a sorted array.
        
    SCHEME:
    quickSort(array, leftmostIndex, rightmostIndex)
      if (leftmostIndex < rightmostIndex)
        pivotIndex <- partition(array,leftmostIndex, rightmostIndex)
        quickSort(array, leftmostIndex, pivotIndex - 1)
        quickSort(array, pivotIndex, rightmostIndex)

    partition(array, leftmostIndex, rightmostIndex)
      set rightmostIndex as pivotIndex
      storeIndex <- leftmostIndex - 1
      for i <- leftmostIndex + 1 to rightmostIndex
        if element[i] < pivotElement
          swap element[i] and element[storeIndex]
          storeIndex++
        swap pivotElement and element[storeIndex+1]
    return storeIndex + 1

    APPLICATIONS:
    Quicksort algorithm is used when:
       the programming language is good for recursion
       time complexity matters
       space complexity matters"""
    
    #Function to divide the array in two parts (left the lower numbers than the 
    #last one (pivot); right (higher ones)) placing pivot in the right index
    def partition(ARRAY, lowestI, highestI, draw_array):
        pivot = ARRAY[highestI]
        i = lowestI - 1
        
        for j in range(lowestI, highestI):
            if ARRAY[j] <= pivot:
                i += 1
                temp = ARRAY[i]
                ARRAY[i] = ARRAY[j]
                ARRAY[j] = temp
                if draw_array:
                    draw_array()
        
        temp = ARRAY[i+1]
        ARRAY[i+1] = ARRAY[highestI]
        ARRAY[highestI] = temp
        if draw_array:
                draw_array()
        
        return i+1
    
    #If we do not introduce a upper limit for index we set it as the lenght
    if highestIndex == False:
        highestIndex = len(array) - 1
        
    #Similar to mergeSort but first divide and then sort
    if lowestIndex < highestIndex:
        pivotIndex = partition(array, lowestIndex, highestIndex, draw_array)
        array = QuickSort(array, lowestIndex, pivotIndex - 1, draw_array)
        array = QuickSort(array, pivotIndex + 1, highestIndex, draw_array)

    return array


def CountingSort(array, draw_array = None):
    """              Counting sort: O(n + k); Stable
    
    Counting sort is a sorting algorithm that sorts the elements of an array 
    by counting the number of occurrences of each unique element in the array.
    The count is stored in an auxiliary array and the sorting is done by 
    mapping the count as an index of the auxiliary array.
    
    SCHEME:
    countingSort(array, size)
      max <- find largest element in array
      initialize count array with all zeros
      for j <- 0 to size
        find the total count of each unique element and 
        store the count at jth index in count array
      for i <- 1 to max
        find the cumulative sum and store it in count array itself
      for j <- size down to 1
        restore the elements to array
        decrease count of each element restored by 1
        
    APPLICATIONS:
    Counting sort is used when: (Only positive integers)
        there are smaller integers with multiple counts.
        linear complexity is the need."""
        
    size = len(array)
    
    output = [0]*size
    #Finding maximum and checking all numbers are positive integers
    max_i = 0
    for i in range(size):
        if array[i]>array[max_i]:
            max_i = i
        if array[i] < 0 or type(array[i]) != int:
            print("[ERROR] Counting Sort only works for positive integers")
            return array
    maximum = array[max_i]
    #Initialize counting array
    count = [0]*(maximum + 1)
    #Fill count array
    for i in range(size):
        count[array[i]] += 1
    for i in range(1, maximum+1):
        count[i] += count[i-1]
    i = size - 1
    while i >= 0:
        output[count[array[i]] - 1] = array[i]
        count[array[i]] -= 1
        i -= 1
        if draw_array:
                draw_array()
    for i in range(size):
        array[i] = output[i]
        if draw_array:
                draw_array()
    
    return array
        
    
def RadixSort(array, draw_array = None):
    """              Radix sort: O(n + k); Stable
    
    Radix sort is a sorting algorithm that sorts the elements by first grouping
    the individual digits of the same place value. Then, sort the elements 
    according to their increasing/decreasing order.
    Suppose, we have an array of 8 elements. First, we will sort elements based
    on the value of the unit place. Then, we will sort elements based on the 
    value of the tenth place. This process goes on until the last significant 
    place.
    
    SCHEME:
    radixSort(array)
       d <- maximum number of digits in the largest element
       create d buckets of size 0-9
       for i <- 0 to d
           sort the elements according to ith place digits using countingSort    
    
    APPLICATIONS:
    Radix sort is implemented in: (Only positive integers (uses counting sort))
       DC3 algorithm (Kärkkäinen-Sanders-Burkhardt) while making a suffix array
       places where there are numbers in large ranges."""
    def countingSort_place(array, place, draw_array):
        size = len(array)
        
        output = [0]*size
        #Finding maximum
        maximum = max(array)
        #Initialize counting array
        count = [0]*(maximum + 1)
        #Fill count array
        for i in range(size):
            if array[i] < 0 or type(array[i]) != int:
                print("[ERROR] Counting Sort only works for positive integers")
                return array
            index = array[i]//place
            count[index % 10] += 1
        for i in range(1, maximum+1):
            count[i] += count[i-1]
        i = size - 1
        while i >= 0:
            index = array[i]//place
            output[count[index % 10] - 1] = array[i]
            count[index % 10] -= 1
            i -= 1
        for i in range(size):
            array[i] = output[i]
            if draw_array:
                draw_array()
        
        return array
    
    max_element = max(array)
    place = 1
    while max_element // place > 0:
        array = countingSort_place(array, place, draw_array)
        place *= 10
    return array

def BucketSort(array, draw_array = None):
    """              Bucket sort: O(n); Stable
    
    Bucket Sort is a sorting algorithm that divides the unsorted array elements
    into several groups called buckets. Each bucket is then sorted by using any
    of the suitable sorting algorithms or recursively applying the same bucket 
    algorithm.Finally, the sorted buckets are combined to form a final sorted 
    array.
    
    
    SCHEME:
    bucketSort()
        create N buckets each of which can hold a range of values
        for all the buckets
            initialize each bucket with 0 values
        for all the buckets
            put elements into buckets matching the range
        for all the buckets 
            sort elements in each bucket
        gather elements from each bucket
    end bucketSort    
        
    APPLICATIONS:
    Bucket sort is used when:
            input is uniformly distributed over a range.
            there are floating point values"""
    size = len(array)
    N = size #Number of buckets
    buckets = []
    [buckets.append([]) for _ in range(N)] #Create N buckets
    maximum = max(array)
    minimum = min(array)
    d = (maximum - minimum)/N #Range of each bucket
    #Fill each bucket with the correspondent element
    for element in array:
        index_bucket = int((element-minimum)/d)
        if element == maximum:#The maximum is out of range, so place it at last bucket
            index_bucket -= 1
        buckets[index_bucket].append(element)
    #Sort each bucket
    for i in range(N):
        buckets[i] = sorted(buckets[i])
    #Combine all the buckets
    index = 0
    for i in range(size):
        for j in range(len(buckets[i])):
            array[index] = buckets[i][j]
            index += 1
            if draw_array:
                draw_array()
        
    return array


def HeapSort(array, draw_array = None):
    """              Heap sort: O(n · log n); NON-Stable
    
    Heap Sort is a popular and efficient sorting algorithm in computer progra-
    mming. Learning how to write the heap sort algorithm requires knowledge of
    two types of data structures - arrays and trees.
    
    SCHEME:
    1. Since the tree satisfies Max-Heap property, then the largest item is 
    stored at the root node.
    2. Swap: Remove the root element and put at the end of the array (nth 
    position) Put the last item of the tree (heap) at the vacant place.
    3. Remove: Reduce the size of the heap by 1.
    4. Heapify: Heapify the root element again so that we have the highest
    element at root.
    5. The process is repeated until all the items of the list are sorted.
        
    
    APPLICATIONS:
    Although Heap Sort has O(n log n) time complexity even for the worst case, 
    it doesn't have more applications ( compared to other sorting algorithms
    like Quick Sort, Merge Sort )."""
    def heapify(array, n, i, draw_array):
        #Find the largest among children
        largest = i
        left_child = 2*i+1
        right_child = 2*i+2
        if left_child < n and array[largest] < array[left_child]:
            largest = left_child
        if right_child < n and array[largest] < array[right_child]:
            largest = right_child
        #If root is not largest, swap and end the heapify
        if largest != i:
            array[i], array[largest] = array[largest], array[i]
            if draw_array:
                draw_array()
            array = heapify(array, n, largest, draw_array) #Heapify the node swaped
            
        return array
    
    
    n = len(array)
    #Make the tree a max heap
    for i in range(n//2, -1, -1):
        array = heapify(array, n, i, draw_array)
    
    for i in range(n-1, 0, -1):
        #Swap root and last index
        array[0], array[i] = array[i], array[0]
        if draw_array:
            draw_array()
        #Heapify the root element
        array = heapify(array, i, 0, draw_array)
        
    return array


def ShellSort(array, draw_array = None):
    """              Shell sort: O(n · log n); NON-Stable
    
    Shell sort is a generalized version of the insertion sort algorithm. It 
    first sorts elements that are far apart from each other and successively 
    reduces the interval between the elements to be sorted.The interval between
    the elements is reduced based on the sequence used.
    
    SCHEME:
    shellSort(array, size)
        for interval i <- size/2n down to 1
            for each interval "i" in array
                sort all the elements at interval "i"
    end shellSort        
        
    APPLICATIONS:
    Shell sort is used when:
        calling a stack is overhead. uClibc library uses this sort.
        recursion exceeds a limit. bzip2 compressor uses it.
        Insertion sort does not perform well when the close elements are far 
        apart. Shell sort helps in reducing the distance between the close 
        elements. Thus, there will be less number of swappings to be performed.    
    """
    n = len(array)
    interval = n//2
    while interval > 0:
        for i in range(interval, n):
            temp = array[i]
            j = i
            while j >= interval and array[j-interval] > temp:
                array[j] = array[j-interval]
                j -= interval
                if draw_array:
                    draw_array()
            array[j] = temp
            if draw_array:
                draw_array()
        interval //= 2
    return array