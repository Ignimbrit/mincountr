#' Load image from file or URL
#' 
#' @description simple wrapper around \code{\link[imager]{load.image}} included
#' here for convenience and consitency reasons only.
#' 
#' @param file path to file of URL
#' 
#' @return an object of class \code{\link[imager]{cimg}}
#' 
#' @examples 
#' mcr_load_image(system.file("extdata", "testim.png", package = "mincountr"))
#' 
#' @export
#' 
mcr_load_image <- function(file){
  imager::load.image(file)
}

#' Plot an images brightness density distribution
#' 
#' @description Inspecting the density distribution of an electron 
#' microscope image of minerals effectively allows to identify different phases
#' as their distinct brightness is a key giveaway of their different chemical
#' composition (and conductivity).
#' 
#' @param x an image loaded with \code{\link{mcr_load_image}}
#' 
#' @examples 
#' myimage <- mcr_load_image(system.file("extdata", "testim.png", package = "mincountr"))
#' mcr_inspect_phases(myimage)
#' 
#' @export

mcr_inspect_phases <- function(x){
  imgr <- imager::grayscale(x) #Making sure image is in greyscale
  imtb <- as.data.frame(imgr)
  imtidy <- tibble::as.tibble(imtb)
  ggplot2::ggplot(data = imtidy, ggplot2::aes(value)) +
    ggplot2::geom_density(kernel="gaussian") +
    ggplot2::scale_x_continuous(breaks=c(seq(0,1,0.1)))
}

#' Plot group assignment of a phase in an image
#' 
#' @description Plots a false color image using groups of brightness that can 
#' either be determined from a graphical inspection of the result of
#' \code{\link{mcr_inspect_phases}} or automaticalle determined by
#' \code{\link{mcr_herd_minerals}}
#' 
#' @param x an image loaded with \code{\link{mcr_load_image}}
#' @param lhs a vector with the left-hand-side position of peaks observed in \code{\link{mcr_inspect_phases}}
#' @param rhs a vector with the right-hand-side position of peaks observed in \code{\link{mcr_inspect_phases}}
#' 
#' @examples 
#' myimage <- mcr_load_image(system.file("extdata", "testim.png", package = "mincountr"))
#' 
#' # (Semi-)Automatic approach:
#' mypeaks <- mcr_autoconstrain(myimage)
#' mcr_inspect_assignment(myimage, mypeaks$x1, mypeaks$x2)
#' 
#' # Manual assignments of brightnessgrouplimits:
#' mcr_inspect_assignment(
#' myimage,
#' lhs = c(0, 0.3, 0.5, 0.92),
#' rhs = c(0.05, 0.45, 0.65, 1)
#' )
#' 
#' @export

mcr_inspect_assignment <-  function(x, lhs, rhs) {
  
  #Some input testing
  if(!all(is.numeric(c(lhs, rhs)))){
    stop("both 'lhs' and 'rhs' must be vectors of type 'numeric'")
  }
  
  if(length(lhs) != length(rhs)){
    stop("numeric vectors 'lhs' and 'rhs' must be of the same length")
  }
  
  if(any(c(lhs, rhs) < 0)){
    stop("elements in numeric vectors 'lhs' ans 'rhs' must be >= 0 & <= 1")
  }
  
  if(any(c(lhs, rhs) > 1)){
    stop("elements in numeric vectors 'lhs' ans 'rhs' must be >= 0 & <= 1")
  }
  
  imgr <- imager::grayscale(x) #Making sure image is in greyscale
  imtb <- as.data.frame(imgr)
  imtidy <- tibble::as.tibble(imtb)
  
  minbins <- as.vector(rbind(lhs, rhs)) 
  
  # assining image pixels to phases, respectively their phases 
  minlev<-.bincode(imtidy$value,minbins,TRUE,TRUE) 
  lvldimg<-cbind(imtidy,minlev)
  
  ggplot2::ggplot(
    data = lvldimg,
    ggplot2::aes(
      x = x, y = y, fill = as.factor(minlev)
    )
    ) +
    ggplot2::geom_raster() +
    ggplot2::labs(fill = "Phase_ID") +
    paletteer::scale_fill_paletteer_d("ggsci::default_igv")
}

#' Automatic peak detection
#' 
#' @description This function automatically determines the position of peaks 
#' and their half-height-width in brightness density distribution of an image.
#' This is key input information necessary to assign brightness-niveaus to
#' certain (mineral) phases an can be used in \code{\link{mcr_herd_minerals}}
#' 
#' @param x an image loaded with \code{\link{mcr_load_image}}
#' 
#' @return A \code{\link[tibble]{tibble}} containing the positions of peak maxima and half-height-width both to the left and to the right of the peak
#' 
#' @examples 
#' myimage <- mcr_load_image(system.file("extdata", "testim.png", package = "mincountr"))
#' mcr_autoconstrain(myimage)
#'  
#' @export
#' 

mcr_autoconstrain <- function(x){
  imgr <- imager::grayscale(x) #Making sure image is in greyscale
  imtb <- as.data.frame(imgr)
  imtidy <- tibble::as.tibble(imtb)
  
  img_density <- stats::density(imtidy$value)
  
  # identifying turning points in the density line
  extreme_points <- tibble::tibble(
    x = img_density$x, y = img_density$y, 
    extr = c(0, diff(sign(diff(img_density$y))), 0)
  )
  
  peaks <- extreme_points %>% 
    dplyr::select(x, extr) %>% 
    dplyr::filter(extr == -2) %>% 
    .$x
  
  peaks[peaks<0] <- 0
  peaks[peaks>1] <- 1
  
  valleys<- extreme_points %>% 
    dplyr::select(x, extr) %>% 
    dplyr::filter(extr == 2) %>% 
    .$x
  
  valleys[valleys<0] <- 0
  valleys[valleys>1] <- 1
  
  # processing each peak individually
  catcher <- vector("list", length = length(peaks))
  for(i in seq_along(catcher)){
    
    # cut out the single peak
    if(peaks[i] == 0){
      cutoff_left <- 0
      cutoff_right <- valleys[valleys > peaks[i]] %>% min()
      if(is.infinite(cutoff_right)){cutoff_right <- 1}
    } else if (peaks[i] == 1){
      cutoff_right <- 1
      cutoff_left <- valleys[valleys < peaks[i]] %>% max()
      if(is.infinite(cutoff_left)){cutoff_left <- 0}
    } else {
      cutoff_left <- valleys[valleys < peaks[i]] %>% max()
      if(is.infinite(cutoff_left)){cutoff_left <- 0}
      cutoff_right <- valleys[valleys > peaks[i]] %>% min()
      if(is.infinite(cutoff_right)){cutoff_right <- 1}
    }
    df <- extreme_points %>% 
      dplyr::filter(x >= cutoff_left) %>% 
      dplyr::filter(x <= cutoff_right) %>% 
      dplyr::mutate(Peak_id = i)
    
    # calculate half-height-width
    
    # left "border" of the actual peak
    if(cutoff_left == 0){
      x1 <- 0
    } else {
        x1 <- df$x[df$x < peaks[i]][which.min(abs(df$y[df$x < peaks[i]]-max(df$y)/2))]
    }
    
    # right "border" of the actual peak
    if(cutoff_right == 1){
      x2 <- 1
    } else {
      x2 <- df$x[df$x > peaks[i]][which.min(abs(df$y[df$x > peaks[i]]-max(df$y)/2))]
    }
    
    catcher[[i]] <- tibble::tibble(
      x1 = x1,
      peakpos = peaks[i],
      x2 = x2,
      ID = i
    )
    
  }
  do.call("rbind", catcher)
}


#' Calculate phase area share
#' 
#' @description calculate the percentage share of mineral phases of the area of
#' an image.
#' 
#' @param x an image loaded with \code{\link{mcr_load_image}}
#' @param lhs a vector with the left-hand-side position of peaks observed in \code{\link{mcr_inspect_phases}}
#' @param rhs a vector with the right-hand-side position of peaks observed in \code{\link{mcr_inspect_phases}}
#'
#' @return A \code{\link[tibble]{tibble}} containing the number of pixels and their relative share of total image area per phase
#'   
#' @examples
#' myimage <- mcr_load_image(system.file("extdata", "testim.png", package = "mincountr"))
#' mypeaks <- mcr_autoconstrain(myimage)
#' mcr_herd_minerals(myimage, mypeaks$x1, mypeaks$x2)
#' 
#' @export
#' 

mcr_herd_minerals <- function(x, lhs, rhs){
  
  #Some input testing
  if(!all(is.numeric(c(lhs, rhs)))){
    stop("both 'lhs' and 'rhs' must be vectors of type 'numeric'")
  }
  
  if(length(lhs) != length(rhs)){
    stop("numeric vectors 'lhs' and 'rhs' must be of the same length")
  }
  
  if(any(c(lhs, rhs) < 0)){
    stop("elements in numeric vectors 'lhs' ans 'rhs' must be >= 0 & <= 1")
  }
  
  if(any(c(lhs, rhs) > 1)){
    stop("elements in numeric vectors 'lhs' ans 'rhs' must be >= 0 & <= 1")
  }
  
  imgr <- imager::grayscale(x) #Making sure image is in greyscale
  imtb <- as.data.frame(imgr)
  imtidy <- tibble::as.tibble(imtb)
  
  
  minbins <- as.vector(rbind(lhs, rhs)) 
  
  # assining image pixels to phases, respectively their phases 
  minlev<-.bincode(imtidy$value,minbins,TRUE,TRUE) 
  lvldimg<-cbind(imtidy,minlev)
  
  sumup <- lvldimg %>%
    dplyr::group_by(minlev) %>%
    dplyr::summarize(count=dplyr::n()) %>% 
    dplyr::ungroup() %>% 
    dplyr::rename(Phase_ID = minlev, pixels = count) %>% 
    dplyr::mutate(
      `proportion_percentage` = (pixels/sum(pixels))*100
    )
  
  return(sumup)
}

