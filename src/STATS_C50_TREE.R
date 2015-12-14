#/***********************************************************************
# * Licensed Materials - Property of IBM 
# *
# * IBM SPSS Products: Statistics Common
# *
# * (C) Copyright IBM Corp. 2015
# *
# * US Government Users Restricted Rights - Use, duplication or disclosure
# * restricted by GSA ADP Schedule Contract with IBM Corp. 
# ************************************************************************/

# author__ = "SPSS, JKP"
# version__ = "1.0.0"

# History
# 22-jul-2015 Original Version
# 15-oct-2015 Treat bands == 1 as zero bands



### MAIN ROUTINE ###
doc50 = function(dep=NULL, indep=NULL, rules=FALSE, missing="include", trials=1,
        usermissingvalid=FALSE,
        subset=TRUE, bands=0, winnow=FALSE, CF=25, minCases=2, costs=NULL,
        samplef=0, seed=NULL, rmode="estimate", useworkspace=FALSE,
        workspaceinfile=NULL, id=NULL, varimp=TRUE, predtype="class",
        factorconversion="levels", outtree="TRUE", outcost=TRUE, influence=NULL,
        workspaceaction="clear", workspacefile=NULL, dataset=NULL) {
    #estimate or predict C5.0 tree
    
    setuplocalization("STATS_C50_TREE")
    
    # A warnings proc name is associated with the regular output
    # (and the same omsid), because warnings/errors may appear in
    # a separate procedure block following the regular output
    procname=gtxt("C5.0 Tree Model")
    warningsprocname = gtxt("C5.0 Tree Model: Warnings")
    omsid="STATSC50"
    warns = Warn(procname=warningsprocname,omsid=omsid)

    tryCatch(library(C50), error=function(e){
        warns$warn(gtxtf("The R %s package is required but could not be loaded.", "C50"),dostop=TRUE)
        }
    )

    if (is.null(seed)) {
        seed = sample.int(4096, size = 1) - 1L
    }
    if (rmode == "estimate" && (is.null(dep) || is.null(indep))) {
        warns$warn(gtxt("Dependent and independent variables must be specified in estimation mode"),
            dostop=TRUE)
    }
    if (rmode == "estimate" && (useworkspace || !is.null(workspaceinfile))) {
        warns$warn(gtxt("Input workspace parameters do not apply to estimation mode"),
            dostop=TRUE)
    }
    if (rmode == "predict" && (!is.null(dep) || !is.null(indep))) {
        warns$warn(gtxt("Dependent and independent variables cannot be specified in prediction mode"),
            dostop=TRUE)
    }
    if (rmode == "predict" && (!useworkspace && is.null(workspaceinfile))) {
        warns$warn(gtxt("A workspace file or retained workspace must be used in prediction mode"),
            dostop=TRUE)
    }
#     if (rmode == "predict" && (workspaceaction== "retain" || !is.null(workspacefile))) {
#         warns$warn(gtxt("The model workspace cannot be retained or saved in prediction mode"),
#             dostop=TRUE)
#     }
    if (rmode =="predict" && is.null(dataset)) {
        warns$warn(gtxt("A dataset name must be supplied if doing predictions"),
            dostop=TRUE)
    }

    checkdataset(dataset, warns)
    
    na.action = ifelse(missing=="omit", na.omit, na.pass)
    # a zero in samplef measns all cases used for training.
    trainprop = samplef
    if (samplef == 1.) {
        samplef = 0
    }
    if (!is.null(spssdictionary.GetWeightVariable())) {
        warns$warn(
        gtxt("The dataset is weighted, but case weights are not used in this procedure except for screening out cases with a nonpositive weight"),
            dostop=FALSE)
    }
    allvars = c(indep, influence, dep)


    if (rmode == "estimate") {
        spssdict = spssdictionary.GetDictionaryFromSPSS(allvars)
        #ylabel = spssdict[length(allvars)- (!is.null(id))]["varLabel",]
        ylabel = spssdict[length(allvars)]["varLabel",]
        if (ylabel=="") {
            ylabel = dep
        }
        # C50imp fails if label contains any regex metacharacters
        # want to just escape them, but backreference not working here
        ylabel = gsub("([][\\()|^$+.])", " ", ylabel, fixed=FALSE)
        allargsest= as.list(environment())
        allargsest[["estdate"]] = date()
        allargsest$packageVersion = packageVersion("C50")
        allargspred = NULL
        # 1 band is equivalent to 0 bands, but the package does not accept that value
        if (bands == 1) {
            bands = 0
        }
    }
    if (rmode == "predict") {
        if (!is.null(workspaceinfile)) {
            load(workspaceinfile)
        }
        # class of a nonexistent object will be NULL
        if (!exists("allargsest") || class(allargsest[["res"]]) != "C5.0") {
            warns$warn(gtxt("The workspace file or retained workspace does not contain a C5.0 model"),
                dostop=TRUE)
        }
        allargspred = as.list(environment())
        dta = spssdata.GetDataFromSPSS(allargsest$indep, factorMode=factorconversion,
            keepUserMissing=usermissingvalid, missingValueToNA=TRUE, row.label=id)

        # check predictor factor status vs estimation time
        predfactors = sapply(dta[1:length(allargsest$indep)], is.factor)
        predfactors = names(predfactors[predfactors])
        factordiffs = setdiff(predfactors, allargsest$factorlist)
        if (length(factordiffs) > 0) {
            warns$warn(gtxtf("Mesurement levels in the prediction data are inconsistent with the estimated model for variables %s",
                paste(factordiffs, collapse=", ")), dostop=TRUE)
        }
    }
    if (rmode == "estimate") {
        dta = spssdata.GetDataFromSPSS(allvars, factorMode=factorconversion, 
            keepUserMissing=usermissingvalid, missingValueToNA=TRUE, row.label=id)
        if (!is.factor(dta[[dep]])) {
            warns$warn("The dependent variable must be categorical", dostop=TRUE)
        }
        if (!is.null(influence) && is.factor(dta[[influence]])) {
            warns$warn(gtxt("The influence variable cannot be categorical"),
                dostop=TRUE)
        }
        # record which predictors are factors
        estfactors = sapply(dta[1:length(indep)], is.factor)
        allargsest$factorlist = names(estfactors[estfactors])
        costmatrix = buildcost(dta[[dep]], costs, warns)
        
        allargsest[["costmatrix"]] = costmatrix
        control = C5.0Control(
            subset=subset,
            bands=bands,
            winnow=winnow,
            CF=CF/100.,
            minCases=minCases,
            sample=samplef,
            seed=seed,
            label = ylabel
        )

        if (is.null(influence)) {
            weights = NULL
        } else {
            weights = dta[[influence]]
        }

        res = tryCatch(C5.0(x=dta[1:length(indep)], 
                y=dta[[dep]], 
                trials=trials,
                weights=weights,
                rules=rules,
                costs=costmatrix,
                na.action=na.action,
                control=control),
                error=function(e) warns$warn(e, dostop=TRUE)
        )
        allargsest[["res"]] = res

    }

    if (!is.null(dataset)) {
        savepred(rmode, allargsest, allargspred, dta, warns)
    }
    ##print(summary(res))

    displayresults(rmode, allargsest, allargspred, dta, warns)
   
    if (workspaceaction == "retain" && rmode=="estimate") {
        assign("allargsest", allargsest, envir=.GlobalEnv)
    }
    if (!is.null(workspacefile)) {
        save(allargsest, file=workspacefile)
    }
    warns$display()
    if (workspaceaction == "clear") {
        rm(list=ls())
        rm(list=ls(envir=.GlobalEnv), envir=.GlobalEnv)
    }
}


buildcost = function(y, costs, warns) {
    # build a cost matrix if any costs were specified
    # y is the dependent variable data
    # costs is the cost specification written as a set of triples
    # costs can have optional = between coordinates and values
    # This is just syntactic suger and is ignored
    #  rowcat, colcat, costn
    # diagonals are ignored
    
    if (is.null(costs)) {
        return(NULL)
    }
    ylevels = levels(y)
    nlevels = length(ylevels)
    costs = costs[costs != "="]
    lencosts = length(costs)

    if (lencosts %% 3 != 0) {
        warns$warn(gtxt("The cost matrix must have a row, column, value specification for each entry"),
            dostop=TRUE)
    }
    costm = matrix(data=rep(NA, nlevels*nlevels), nrow=nlevels)
    rownames(costm) = ylevels
    colnames(costm) = ylevels
    for (i in ylevels) {
        costm[[i,i]] = 0.
    }
    tryCatch(
        {
        for (i in seq(1, lencosts, 3)) {
            costm[[costs[[i]], costs[[i+1]]]] = as.double(costs[i+2])
        }
        },
        error = function(e) warns$warn(gtxtf("A cost matrix index or value is invalid: %s, %s",
                    costs[i], costs[i+1]), dostop=TRUE)
    )
    if (any(is.na(costm))) {
        warns$warn(gtxt("One or more cost matrix values were not specified"), dostop=TRUE)
    }
    return(costm)
}


savepred = function(rmode, allargsest, allargspred, dta, warns) {
    # create a dataset of predicted values
    
    if (rmode == "estimate") {
        dataset = allargsest$dataset
        trials = allargsest$trials
        na.action = allargsest$na.action
        ptype = allargsest$predtype
        idname = allargsest$id
    } else {
        dataset = allargspred$dataset
        trials = allargspred$trials
        na.action = allargspred$na.action
        ptype = allargspred$predtype
        idname = allargspred$id
    }
    if (is.null(dataset)) {
        return()
    }

    if (!is.null(allargsest$costmatrix) && ptype=="prob") {
        warns$warn(gtxt("Probabilities are not available with a cost matrix.  Switching to class"),
            dostop=FALSE)
        ptype="class"
    }
    idlabel = ifelse(is.null(idname), "ID", idname)
    # dta will have the ID variable values as the row names if given
    pred = tryCatch(
            data.frame(predict(allargsest$res, type=ptype, 
                newdata=dta, na.action=na.action, trials=trials)),
            error=function(e) {
                warns$warn(e, dostop=TRUE)
            }
            )

    # create dictionary
    # length for ID values allowing for Unicode expansion
    idlen = 3 * max(sapply(row.names(pred), length))
    dict = list()
    dict[[1]] = c("ID", idlabel, idlen, paste("A", idlen, sep=""), "nominal")
    if (ptype == "class") {
        if (is.character(pred[1,1])) {
            vlen = 3 * max(sapply(pred[[1]], length))
            dict[[2]] = c(
                "PredClass",
                gtxt("Predicted Class"),
                vlen,
                paste("A", vlen, sep=""),
                "nominal"
            )
        } else {
            dict[[2]] = c(
                "PredClass", 
                gtxt("Predicted Class"), 
                0, 
                "F8.0",
                "nominal"
            )
        }
    } else {  #probabilities
        thenames = names(pred)
        for (i in 1:ncol(pred)) {
            dict[[i+1]] = c(
                paste("P", i, sep=""),
                thenames[[i]], 
                0,
                "F8.4", 
                "scale"
            )
        }
    }
    dict = spssdictionary.CreateSPSSDictionary(dict)
    spssdictionary.SetDictionaryToSPSS(dataset, dict)
    spssdata.SetDataToSPSS(dataset, data.frame(row.names(pred), pred))
    spssdictionary.EndDataStep()
}


displayresults = function(rmode, allargsest, allargspred, dta, warns) {
    # Display output
    
    StartProcedure(allargsest[["procname"]], allargsest[["omsid"]])
    summarylabels=list(
        gtxt("Dependent Variable"),
        gtxt("Independent Variables"),
        gtxt("ID Variable"),
        gtxt("Model Source"),
        gtxt("Estimation Date"),
        gtxt("Missing Values"),
        gtxt("User Missing Value Treatment"),
        gtxt("Factor Conversion"),
        gtxt("Influence"),
        gtxt("Boosts Requested"),
        gtxt("Boosts Actual"),
        gtxt("Number of Samples"),
        gtxt("Number of Predictors"),
        gtxt("Confidence Pct"),
        gtxt("Minimum Size"),
        gtxt("Training Proportion Subsample"),
        gtxt("Random Number Seed"),
        gtxt("Prediction Dataset")
    )
    
    # modespecific parameters
    if (rmode == "estimate") {
        idvar = allargsest$id
        ms = gtxt("--NA--")
        dataset = allargsest$dataset
        outcost = allargsest$outcost
        outtree = allargsest$outtree
    } else {
        idvar = allargspred$id
        ms = ifelse(!is.null(allargspred$workspacefile), allargspred$workspacefile,
            gtxt("--Workspace--"))
        dataset = allargspred$dataset
        outcost = allargspred$outcost
        outtree = allargspred$outtree
    }
    
    # specification summary
    tReq = allargsest$res$trials[[1]]
    tAct = ifelse(length(allargsest$res$trials) > 1, allargsest$res$trials[[2]], 1)
    summaryvalues = list(
        allargsest[["dep"]],
        paste(allargsest$indep, collapse=", "),
        ifelse(is.null(idvar), gtxt("--NA--"), idvar),
        ms,
        allargsest$estdate,
        ifelse(allargsest$missing == "include", gtxt("include"), gtxt("omit")),
        ifelse(allargsest$usermissingvalid, gtxt("Valid"), gtxt("Missing")),
        ifelse(allargsest$factorconversion == "levels", gtxt("levels"), gtxt("labels")),
        ifelse(is.null(allargsest$influence), gtxt("--NA--"), allargsest$influence),
        tReq,
        tAct,
        allargsest$res$dims[[1]],
        allargsest$res$dims[[2]],
        allargsest$CF,
        allargsest$minCases,
        allargsest$trainprop,
        allargsest$seed,
        ifelse(is.null(dataset), gtxt("--NA--"), dataset)
    )

    names(summaryvalues) = summarylabels
    summarydf = data.frame(cbind(summaryvalues))
    colnames(summarydf) = gtxt("Values")

    spsspivottable.Display(
        summarydf, 
        title=gtxt("C5.0 Tree Summary"), 
        templateName="STATSC5.0SUMMARY",
        caption=gtxtf("Results computed by R C5.0 package, version %s", 
            allargsest$packageVersion),
        isSplit=FALSE
    )
    
    # cost matrix
    if (outcost && !is.null(allargsest$costmatrix)) {
        costm = data.frame(allargsest$costmatrix)
        names(costm) = row.names(costm)
        spsspivottable.Display(
            costm, 
            title=gtxt("Cost Matrix"),
            templateName="STATSC5.0COSTS",
            isSplit=FALSE)
    }
    
    # variable importance
    if (rmode == "estimate" && allargsest$varimp) {
            imp = data.frame(C5imp(allargsest$res, metric="usage"))
            imp = data.frame(imp, C5imp(allargsest$res, metric="splits"))
            names(imp) = c(gtxt("Usage"), gtxt("Splits"))             
            spsspivottable.Display(
                imp, 
                title = gtxt("Variable Importance"),
                templateName="STATSC5.0VARIMP"
            )
    }
    
    # tree and summary as text
    if (outtree) {
        results=summary(allargsest$res)
    } else {
        results = allargsest$res
    }

    class(results) = "character"
    spss.TextBlock(gtxt("Estimation Results"), results,
        gtxt("Estimation Results"))
    
    spsspkg.EndProcedure()
}

checkdataset = function(dataset, warns) {
    if (!is.null(dataset)) {
        alldatasets = spssdata.GetDataSetList()
        if ("*" %in% alldatasets) {
            warns$warn(gtxt("The active dataset must have a name in order to use this procedure"),
                       dostop=TRUE)
        }
        if (dataset %in% alldatasets) {
            warns$warn(gtxt("The output dataset name must not already be in use"),
                       dostop=TRUE)
        }
    }
}
Warn = function(procname, omsid) {
    # constructor (sort of) for message management
    lcl = list(
        procname=procname,
        omsid=omsid,
        msglist = list(),  # accumulate messages
        msgnum = 0
    )
    # This line is the key to this approach
    lcl = list2env(lcl) # makes this list into an environment

    lcl$warn = function(msg=NULL, dostop=FALSE, inproc=FALSE) {
        # Accumulate messages and, if dostop or no message, display all
        # messages and end procedure state
        # If dostop, issue a stop.

        if (!is.null(msg)) { # accumulate message
            assign("msgnum", lcl$msgnum + 1, envir=lcl)
            # There seems to be no way to update an object, only replace it
            m = lcl$msglist
            m[[lcl$msgnum]] = msg
            assign("msglist", m, envir=lcl)
        } 

        if (is.null(msg) || dostop) {
            lcl$display(inproc)  # display messages and end procedure state
            if (dostop) {
                stop(gtxt("End of procedure"), call.=FALSE)  # may result in dangling error text
            }
        }
    }
    
    lcl$display = function(inproc=FALSE) {
        # display any accumulated messages as a warnings table or as prints
        # and end procedure state, if any

        if (lcl$msgnum == 0) {   # nothing to display
            if (inproc) {
                spss.EndProcedure()
            }
        } else {
            if (!inproc) {
                procok =tryCatch({
                    StartProcedure(lcl$procname, lcl$omsid)
                    TRUE
                    },
                    error = function(e) {
                        FALSE
                    }
                )
            }
            if (procok) {  # build and display a Warnings table if we can
                table = spss.BasePivotTable("Warnings ","Warnings") # do not translate this
                rowdim = BasePivotTable.Append(table,Dimension.Place.row, 
                    gtxt("Message Number"), hideName = FALSE,hideLabels = FALSE)

                for (i in 1:lcl$msgnum) {
                    rowcategory = spss.CellText.String(as.character(i))
                    BasePivotTable.SetCategories(table,rowdim,rowcategory)
                    BasePivotTable.SetCellValue(table,rowcategory, 
                        spss.CellText.String(lcl$msglist[[i]]))
                }
                spsspkg.EndProcedure()   # implies display
            } else { # can't produce a table
                for (i in 1:lcl$msgnum) {
                    print(lcl$msglist[[i]])
                }
            }
        }
    }
    return(lcl)
}

# localization initialization
setuplocalization = function(domain) {
    # find and bind translation file names
    # domain is the root name of the extension command .R file, e.g., "SPSSINC_BREUSCH_PAGAN"
    # This would be bound to root location/SPSSINC_BREUSCH_PAGAN/lang

    fpath = Find(file.exists, file.path(.libPaths(), paste(domain, ".R", sep="")))
    bindtextdomain(domain, file.path(dirname(fpath), domain, "lang"))
} 
# override for api to account for extra parameter in V19 and beyond
StartProcedure <- function(procname, omsid) {
    if (substr(spsspkg.GetSPSSVersion(),1, 2) >= 19) {
        spsspkg.StartProcedure(procname, omsid)
    }
    else {
        spsspkg.StartProcedure(omsid)
    }
}

gtxt <- function(...) {
    return(gettext(...,domain="STATS_C50_TREE"))
}

gtxtf <- function(...) {
    return(gettextf(...,domain="STATS_C50_TREE"))
}


Run = function(args) {
    #Execute the STATS C50 TREE command

    cmdname = args[[1]]
    args = args[[2]]
    oobj = spsspkg.Syntax(list(
        spsspkg.Template("MODE", subc="", ktype="str", var="rmode",
            vallist=list("estimate", "predict")),
        spsspkg.Template("DEPENDENT", subc="",  ktype="existingvarlist", var="dep"),
        spsspkg.Template("INDEPENDENT", subc="", ktype="existingvarlist", var="indep", 
            islist=TRUE),
        spsspkg.Template("ID", subc="", ktype="existingvarlist", var="id"),
        spsspkg.Template("INFLUENCE", subc="", ktype="existingvarlist", var="influence"),
        spsspkg.Template("COSTS", subc="", ktype="literal", var="costs", islist=TRUE),
        spsspkg.Template("USEWORKSPACE", subc="", ktype="bool", var="useworkspace"),
        spsspkg.Template("WORKSPACEFILE", subc="", ktype="literal", var="workspaceinfile"),
        
        spsspkg.Template("FACTORCONVERSION", subc="OPTIONS", ktype="str", var="factorconversion",
            vallist=list("labels", "levels")),
        spsspkg.Template("MISSING", subc="OPTIONS", ktype="str", var="missing", 
            vallist=list("omit", "include")),
        spsspkg.Template("USERMISSINGVALID", subc="OPTIONS", ktype="bool", var="usermissingvalid"),
        spsspkg.Template("BOOSTS", subc="OPTIONS", ktype="int", var="trials",
            vallist=list(1)),
        spsspkg.Template("GROUPING", subc="OPTIONS", ktype="bool", var="subset"),
        spsspkg.Template("FEATURESELECTION", subc="OPTIONS", ktype="bool", var="winnow"),
        spsspkg.Template("CONFIDENCE", subc="OPTIONS", ktype="float", var="CF",
            vallist=list(.01, 99.9999)),
        spsspkg.Template("MINSIZE", subc="OPTIONS", ktype="int", var="minCases"),
        spsspkg.Template("TRAINPROP", subc="OPTIONS", ktype="float", var="samplef",
            vallist=list(0, 1)),
        spsspkg.Template("RNSEED", subc="OPTIONS", ktype="int", var="seed"),
        
        spsspkg.Template("TREE", subc="PRINT", ktype="bool", var="outtree"),
        spsspkg.Template("RULES", subc="PRINT", ktype="bool", 
            var="rules"),
        spsspkg.Template("RULEBANDS", subc="PRINT", ktype="int", var="bands",
            vallist=list(0,1000)),
        spsspkg.Template("COSTMATRIX", subc="PRINT", ktype="bool", var="outcost"),
        spsspkg.Template("VARIMPORTANCE", subc="PRINT", ktype="bool", var="varimp"),
        
        spsspkg.Template("WORKSPACEACTION", subc="SAVE", ktype="str", var="workspaceaction",
            vallist=list("retain", "clear")),
        spsspkg.Template("WORKSPACEOUTFILE", subc="SAVE", ktype="literal", var="workspacefile"),
        spsspkg.Template("DATASET", subc="SAVE", ktype="varname", var="dataset"),
        spsspkg.Template("PREDTYPE", subc="SAVE", ktype="str", var="predtype",
            vallist=list("class", "prob"))
    ))

    # A HELP subcommand overrides all else
    if ("HELP" %in% attr(args,"names")) {
        #writeLines(helptext)
        helper(cmdname)
    }
    else {
        res <- spsspkg.processcmd(oobj, args, "doc50")
    }
}

helper = function(cmdname) {
    # find the html help file and display in the default browser
    # cmdname may have blanks that need to be converted to _ to match the file
    
    fn = gsub(" ", "_", cmdname, fixed=TRUE)
    thefile = Find(file.exists, file.path(.libPaths(), fn, "markdown.html"))
    if (is.null(thefile)) {
        print("Help file not found")
    } else {
        browseURL(paste("file://", thefile, sep=""))
    }
}
if (exists("spsspkg.helper")) {
assign("helper", spsspkg.helper)
}
