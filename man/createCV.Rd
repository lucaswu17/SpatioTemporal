\name{createCV}
\encoding{latin1}
\Rdversion{1.1}
\alias{createCV}

\title{
  Define Cross-Validation Groups
}
\description{
  Creates a matrix that specifies cross-validation schemes.
}
\usage{
createCV(mesa.data.model, groups = 10, min.dist = 0.1,
         random = FALSE, subset=NA, option="all")
}
\arguments{
  \item{mesa.data.model}{
    Data structure holding observations, and information regarding the
    observation locations. See \code{\link{create.data.model}} and
    \code{\link{mesa.data.model}}.
  }
  \item{groups}{
    Number of cross-validation groups, zero gives leave-one-out
    cross-validation.
  }
  \item{min.dist}{
    Minimum distance between locations for them to end up in separate
    groups. Points closer than \code{min.dist} will be forced into the
    same group. A high value for \code{min.dist} can result in fewer
    cross-validation groups than specified in \code{groups}.
  }
  \item{random}{
    If \code{FALSE} repeated calls to the function will return the same
    grouping, if \code{TRUE} repeated calls will give different
    CV-groupings. Ensures that simulation studies are reproducable.
  }
  \item{subset}{
    A subset of locations for which to define the cross-validation
    setup. Only sites listed in \code{subset} are dropped from one of
    the cross-validation groups; in other words sites \emph{not in}
    \code{subset} are used for estimation and preidiction of \emph{all}
    cross-validation groups.

    This option is \emph{ignored} if \code{option!="all"}.
  }
  \item{option}{
    For internal MESA Air usage, see Details below.
  }
}
\details{
  The number of observations left out of each group is can be rather
  uneven; the main goal of \code{createCV} is to create CV-groups such
  that the groups contain roughly the same \emph{number of locations}
  ignoring the number of observations at each location. If there are
  large differences in the number of  observations at differnt locations
  one could use the \code{subset} option to create different
  CV-groupings for different types of locations. The groups can then be
  combined as \cr \code{I.final = I.1 | I.2 | I.3}.

  If \code{random=FALSE} the function initially sets \cr
  \code{set.seed(0, kind = "Mersenne-Twister")}, \cr
  and resets the random-seed using \code{\link{.Random.seed}} and
  \code{\link{set.seed}} before terminating.

  The \code{option} input determines which sites to include in the
  cross-validation. Possible options are \code{"all"}, \code{"fixed"}, 
  \code{"comco"}, \code{"snapshot"} and \code{"home"}.
  \describe{
    \item{\code{"all"}}{Uses all available sites, possibly subset
      according to \code{subset}. The sites will be grouped with sites
      seperated by less than \code{min.dist} being put in the same
      CV-group.
    }
    \item{\code{"fixed"}}{Uses only sites that have \cr
      \code{mesa.data.model$location$type \%in\%
	c("AQS","FIXED")}. Given the subsettting the sites will be
      grouped as for \code{"all"}.
    }
    \item{\code{"home"}}{Uses only sites that have \cr
      \code{mesa.data.model$location$type \%in\% c("HOME")}. Given the
      subsettting the sites will be grouped as for \code{"all"}.
    }
    \item{\code{"comco"}, \code{"snapshot"}}{Uses only sites that have
      \cr \code{mesa.data.model$location$type \%in\% c("COMCO")}.

      The sites will be grouped  together if they are from the same road
      gradient. The road gradients are grouped by studying the name of
      the sites. With "?" denoting one or more letters and "#" denoting
      one or more digits the names are expected to follow "?-?#?#", for
      random sites, and "?-?#?#?" for the gradients (with all but the
      last letter being the same for the entire gradient).
    }
  }
}
\value{
  Return a (number or observations) - by - (groups) logical
  matrix. Each column defines a cross-validation set with the
  \code{TRUE} values marking the observations to be left out.
}
\author{
  \enc{Johan Lindström}{Johan Lindstrom}
}
\seealso{
  See also \code{\link{estimateCV}}, and \code{\link{predictCV}}. 
  
  For computing CV statistics, see also \code{\link{compute.ltaCV}}, 
  \code{\link{predictNaive}}, and for further illustration see
  \code{\link{plotCV}}, \code{\link{CVresiduals.qqnorm}}, 
  and \code{\link{summaryStatsCV}}.
}
\examples{
##load the data
data(mesa.data.model)

##create a matrix with the CV-schemes
I.cv <- createCV(mesa.data.model, groups=10)

##number of observations in each CV-group
colSums(I.cv)

##Which sites belong to which groups?
ID.cv <- lapply(apply(I.cv,2,list),function(x) 
                unique(mesa.data.model$obs$ID[x[[1]]]))
print(ID.cv)

##Note that the sites with distance 0.084<min.dist 
##are grouped together (in group 10).
mesa.data.model$dist[ID.cv[[10]],ID.cv[[10]]]

##Find out which location belongs to which cv group
I.col <- apply(sapply(ID.cv,function(x) mesa.data.model$location$ID
               \%in\% x), 1, function(x) if(sum(x)==1) which(x) else 0)
names(I.col) <- mesa.data.model$location$ID
print(I.col)

##Plot the locations, colour coded by CV-grouping
plot(mesa.data.model$location$long, mesa.data.model$location$lat,
     pch=23+floor(I.col/max(I.col)+.5), bg=I.col, 
     xlab="Longitude",ylab="Latitude")
}
