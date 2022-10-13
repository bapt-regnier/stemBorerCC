# stemBorerCC
R scripts for the assessment of the impacts of global warming on the development of four stemborers in temperate and tropical climates

Several R packages are required : `targets`, `truncnorm`, `devRate`, `nlstools`, `parallel`

The complete temperature datasets were retrieved from the IPCC WGI Interactive Atlas GitHub repository under CC-BY 4.0 license :

> Iturbide, M., Fernández, J., Gutiérrez, J.M., Bedia, J., Cimadevilla, E., Díez-Sierra, J., Manzanas, R., Casanueva, A., Baño-Medina, J., Milovac, J., Herrera, S., Cofiño, A.S., San Martín, D., García-Díez, M., Hauser, M., Huard, D., Yelekci, Ö. (2021) Repository supporting the implementation of FAIR principles in the IPCC-WG1 Atlas. Zenodo, DOI: 10.5281/zenodo.3691645. Available from: https://github.com/IPCC-WG1/Atlas

Note that at least 32GB of RAM are required to run the scripts with default values. The arguments nInd and tempSeq of the function computeDT() in _targets.R file can be modified to reduce memory usage.

To run the scripts using the `targets` package, first install the packages `truncnorm`, `devRate`, `nlstools`, and `parallel`. Then run `tar_make()`.
