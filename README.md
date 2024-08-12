![example branch parameter](https://github.com/jlab/fold-grammars/actions/workflows/c-cpp.yml/badge.svg)

# fold-grammars
Collection of bgap code for RNA folding

For rudimentary documentation, consult https://bibiserv.cebitec.uni-bielefeld.de/fold-grammars/

# Install derived software
## Ubuntu
Register janssenlab/software as additional Ubuntu Personal Package Archive (PPA):
```
sudo add-apt-repository ppa:janssenlab/software
sudo apt update
```
Install software one or multiple software packages:
  - RNAshapes: `sudo apt-get install rnashapes`
  - RNAalishapes: `sudo apt-get install rnaalishapes`
  - pKiss: `sudo apt-get install pkiss`
  - pAliKiss: `sudo apt-get install palikiss`
  - KnotInFrame: `sudo apt-get install knotinframe`
  - RapidShapes: `sudo apt-get install rapidshapes`
  - aCMs: `sudo apt-get install acms`

## conda
You can alternatively install precompiled binaries for Linux and OSX through (bio)conda:
  - RNAshapes: [`conda install bioconda::rnashapes`](https://anaconda.org/bioconda/rnashapes)
  - RNAalishapes: [`conda install bioconda::rnaalishapes`](https://anaconda.org/bioconda/rnaalishapes)
  - pKiss: [`conda install bioconda::pkiss`](https://anaconda.org/bioconda/pkiss)
  - pAliKiss: [`conda install bioconda::palikiss`](https://anaconda.org/bioconda/palikiss)
  - KnotInFrame: [`conda install bioconda::knotinframe`](https://anaconda.org/bioconda/knotinframe)
  - RapidShapes: [`conda install bioconda::rapidshapes`](https://anaconda.org/bioconda/rapidshapes)
  - aCMs: [`conda install bioconda::acms`](https://anaconda.org/bioconda/acms) (no OSX version available)

(We have not (yet) tested the bioconda linux-aarch64 builds!)
