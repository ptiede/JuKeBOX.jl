# JuKeBOX

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://ptiede.github.io/JuKeBOX.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://ptiede.github.io/JuKeBOX.jl/dev)
[![Build Status](https://github.com/ptiede/JuKeBOX.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/ptiede/JuKeBOX.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/ptiede/JuKeBOX.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/ptiede/JuKeBOX.jl)
[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/invenia/BlueStyle)
[![ColPrac: Contributor's Guide on Collaborative Practices for Community Packages](https://img.shields.io/badge/ColPrac-Contributor's%20Guide-blueviolet)](https://github.com/SciML/ColPrac)


# Installation

Currently JuKeBOX is not a registered package and depend on an additional non-registered package. This makes installation a little more complicated than it needs to be and will be fixed in the future. To install JuKeBOX do the following:

  1. Install a forked version of Elliptic.jl from ()[https://github.com/ptiede/Elliptic.jl]. Specifically we need the extension branch which will make JuKeBOX differentiable. To do this go into the Julia repl press `]` to enter Pkg mode and type
    
    
    add "https://github.com/ptiede/Elliptic.jl#extension"
    
 To add the correct branch
    
    
  2. Add JuKeBOX while in the Pkg mode by typing 
      ```julia
      add "https://github.com/ptiede/JuKeBOX.jl"
      ```
    This will download the main branch which should work

    

