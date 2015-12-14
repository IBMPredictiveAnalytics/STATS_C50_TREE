* Encoding: UTF-8.
file handle data/name="c:/spss23/samples/english/employee data.sav".
get file=data.
dataset name emp.

stats c50 tree dependent=jobcat independent=salary minority gender educ.

stats c50 tree dependent=jobcat independent=salary minority gender educ
costs = 1 2 2.5 1 3 3.0 2 1 1.8 2 3 2.8 3 1 .5 3 2 .75
/options factorconversion=levels
/save workspaceoutfile="c:/temp/treeworkspace.Rdata".

stats c50 tree dependent=jobcat independent=salary minority gender educ
costs = 1 2 2.5 1 3 3.0 2 1 1.8 2 3 2.8 3 1 .5 3 2 .75
/options factorconversion=levels
/print tree=yes costmatrix=yes varimportance=yes
/save workspaceoutfile="c:/temp/treeworkspace.Rdata".

stats c50 tree mode=predict workspacefile="c:/temp/treeworkspace.Rdata"
/save dataset=predict.

* no costs.
stats c50 tree dependent=jobcat independent=salary minority gender educ
/options factorconversion=levels
/print tree=yes costmatrix=yes varimportance=yes
/save workspaceoutfile="c:/temp/treeworkspaceNC.Rdata".

stats c50 tree mode=predict workspacefile="c:/temp/treeworkspaceNC.Rdata"
/save dataset=predict.


* no tree, factor labels.
stats c50 tree dependent=jobcat independent=salary minority gender educ
/options factorconversion=labels
/print tree=no costmatrix=no varimportance=no.

stats c50 tree dependent=gender independent=salary.

* string factor with costs.
stats c50 tree dependent=gender independent=salary
costs= 'f' 'm' 2 'm' 'f' 3.

* negative cost.
stats c50 tree dependent=gender independent=salary
costs= 'f' 'm' -2 'm' 'f' 3.

* more options.
stats c50 tree dependent=jobcat independent=salary minority gender educ
/options factorconversion=labels missing=omit
boosts=5 grouping=yes minsize=10 rnseed=1234 featureselection=yes
confidence=95
/print rules=yes.

compute weight=rv.uniform(0, 2).
variable level weight(nominal).
weight by weight.
weight off.
* should warn.
stats c50 tree dependent=jobcat independent=salary minority gender educ.

stats c50 tree dependent=jobcat independent=salary minority gender educ
/options influence=weight.



