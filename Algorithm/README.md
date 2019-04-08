To run our fancy algorithm
  1. Since we are basing all of our work off of the 64-long array, we should not need to alter our library of valid arrays, but if you want to use a different library
    -> Run the AllArrays.m file. This will create new B and B_bad (M, p, L, totalM matrix of all verified configurations) and sub2s & sub2s_bad (the necessary subarray       2 configuration for those M and L values) matrices that will be saved directly into the folder that you are in (the Matlab path           that you are on when you run the script). These matrices should have the same number of rows!!
    
  2. MainAlgorithm.m if a function file where the input is a logic array (1xn) that represents the on (1) and off (0) sensors in the array. 
    -> Later we will have to change some stuff in here to use real data (maybe add more inputs to this function) and apply filtering. 
    
  3. sensors.m is just a 35-long example logical array that I used to test MainAlgorithm.m  
