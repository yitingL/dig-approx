# Approximate Logic Synthesis for Dot-Inverter Graphs
This project is based on the paper:
- Approximate Logic Synthesis for Dot-Inverter Graphs Using Node Merging-enhanced Genetic Algorithm-based Approach



# To Compile
In this directory,
```
$ make
```

# To Run
In this directory,
```
$ ./main ...
        -g <golden_circuit>
        -o <output_circuit>
        -e <error_rate in %>
        [-t <timeout> = 3600]
        [-c <original_circuit>]
        [-d <ed_mode>]
```

# References
- [ABC: System for Sequential Logic Synthesis and Formal Verification](https://github.com/berkeley-abc/abc)
- https://github.com/lsils/mockturtle