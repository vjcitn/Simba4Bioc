# Simba4Bioc
Elementary Bioconductor interfaces to Pinello Lab [SIMBA graph embeddings](https://www.nature.com/articles/s41592-023-01899-8) for single cell omics

This package provides R functions to build and interact with
PyTorch biggraph artifacts generated using the simba
python module.

The basic vignette (use 'Get Started' tab above) 
shows how to build interactive visualizations
of joint embedding of cell and gene profiles like this one:

<iframe src="reference/figures/litpl.html" width="100%" height="360px"/>

Additional developments will target exploration of multimodal
single cell assays.

(N.B.  There are extraneous "Aa" symbols and mappings in the
plot displayed above.  This is a consequence of the introduction
of extra symbols on the display.  Improved plotly infrastructure
or usage may help remove these meaningless tokens.  Suggestions
are welcome.)

