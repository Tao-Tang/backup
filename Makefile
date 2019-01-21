All:
        icc -std=c++11 -o ECC ECC.cpp -pthread


clean:
        $(RM) $(OUTPUT)
