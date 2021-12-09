System{
  Mixture{
    nMonomer  2
    monomers  A   1.0  
              B   1.0 
    nPolymer  1
    Polymer{
       type    linear
       nBlock  2
       blocks  0    0.3   
               1    0.7  
       phi     1.0
    }
    ds   0.01
  }
  ChiInteraction{
    chi  0   0   0.0
         1   0   20.0
         1   1   0.0
  }

     unitCell    hexagonal   1.6651597308
     mesh        32    32
     groupName   p_6_m_m

  AmIterator{
    maxItr 100
    epsilon 1e-10
    maxHist 15
    isFlexible 1
  }
}