System{
  Mixture{
    nMonomer  2
    monomers[  
              1.0
              1.0
    ]
    nPolymer  1
    Polymer{
       type    linear
       nBlock  2
       blocks[
               0    0.3
               1    0.7
       ]
       phi     1.0
    }
    ds   0.01
  }
  Interaction{
    chi(  
         1   0   20.0
    )
  }
  Domain{
    mesh        32    32
    lattice     hexagonal
    groupName   p_6_m_m
  }
  AmIterator{
    epsilon      1.0e-09
    maxItr       100
    maxHist      30
    verbose      1
    isFlexible   1
  }
}

    lattice     hexagonal
    unitCell    hexagonal   1.6908668724
