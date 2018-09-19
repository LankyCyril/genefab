LIMMA_SCRIPT = """
if (!require("affy")) {{
    source("https://bioconductor.org/biocLite.R")
    biocLite("affy", suppressUpdates=TRUE, suppressAutoUpdate=TRUE, ask=FALSE)
}}
library("affy")

if (!require("limma")) {{
    source("https://bioconductor.org/biocLite.R")
    biocLite("limma", suppressUpdates=TRUE, suppressAutoUpdate=TRUE, ask=FALSE)
}}
library("limma")

affy_set = ReadAffy(filenames=c({cels}))
annotation = data.frame({factor_name}=factor(c({factor_list})))
rownames(annotation) = sampleNames(affy_set)
phenoData(affy_set) = AnnotatedDataFrame(annotation)

expression_set = rma(affy_set)
design = model.matrix(~expression_set${factor_name})
bayes_data = eBayes(lmFit(expression_set, design))
results = topTable(bayes_data, n=Inf)

write.csv(results, file={output}, quote=FALSE)
"""
