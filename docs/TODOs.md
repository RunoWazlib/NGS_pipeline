## TODOs ##
1. Add trimming functions
2. Add the capability to only pass in an alignment and get the plots/analysis out
3. Add secondary structure prediction
    - Global analysis of sequencing files
        - Add file layer / sample treatment to config
    - Mutational frequency plots from SHAPE data & determine modified positions
    - Parameterize modified positions into energy constraints
    - Model secondary structures w RNAStructure `Fold` tool using these constraints
        - I suppose one could do this anyway regardless of SHAPE data or no... cross that bridge later