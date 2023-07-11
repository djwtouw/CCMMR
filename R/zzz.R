.onUnload <- function (libpath)
{
    library.dynam.unload("CCMMR", libpath)
}


.onAttach <- function(libname, pkgname)
{
    package_citation = paste("\nTouw, D.J.W., Groenen, P.J.F., and Terada, Y.",
                             "(2022). Convex Clustering through MM: An",
                             "Efficient Algorithm to Perform Hierarchical",
                             "Clustering. arXiv preprint arXiv:2211.01877.\n")
    message("Thank you for using CCMMR!")
    message("To acknowledge our work, please cite the paper:")
    message(package_citation)
}
