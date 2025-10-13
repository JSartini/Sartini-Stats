
Fourier_bases <- function(Q){
  basis = map(2:Q, function(x){
    flag = x %% 2
    mult = x %/% 2
    if(flag == 0){ # Cosine
      return_func <- function(x){
        sqrt(2)*cos(2*pi*mult*x)
      }
    }
    else{ # Sine
      return_func <- function(x){
        sqrt(2)*sin(2*pi*mult*x)
      }
    }
    return(return_func)
  })
  
  names(basis) = paste0("Func ", 1:(Q-1))
  basis[["Func 0"]] = function(x){return(rep(1, length.out = length(x)))}
  return(basis)
}

Fourier_d2 <- function(Q){
  derivs = map(2:Q, function(x){
    flag = x %% 2
    mult = x %/% 2
    if(flag == 0){ # Cosine
      return_func <- function(x){
        -sqrt(2)*cos(2*pi*mult*x)*(2*pi*mult)^2
      }
    }
    else{ # Sine
      return_func <- function(x){
        -sqrt(2)*sin(2*pi*mult*x)*(2*pi*mult)^2
      }
    }
    return(return_func)
  })
  
  names(derivs) = paste0("Func ", 1:(Q-1))
  derivs[["Func 0"]] = function(x){return(rep(0, length.out = length(x)))}
  
  return(derivs)
}

Legendre_bases <- function(Q){
  basic_funcs = polynomial.functions(legendre.polynomials(Q-1))
  
  scaled_funcs = map(1:length(basic_funcs), function(idx){
    scaling_factor = sqrt(2*(idx-1) + 1)
    out_func <- function(x){
      return(scaling_factor * basic_funcs[[idx]](2*x-1))
    }
  })
  
  return(scaled_funcs)
}

Legendre_d2 <- function(Q){
  basic_derivs = legendre.polynomials(Q-1) %>%
    polynomial.derivatives() %>%
    polynomial.derivatives() %>%
    polynomial.functions()
  
  scaled_derivs = map(1:length(basic_derivs), function(idx){
    scaling_factor = sqrt(2*(idx-1) + 1)
    out_func <- function(x){
      return(scaling_factor * basic_derivs[[idx]](2*x-1) * 4) # second order chain rule
    }
  })
}

Splinet_bases <- function(Q){
  B_f = splinet(knots = seq(0, 1, length.out = Q+2), norm = T)
  basis = map(1:(Q-2), function(q){
    func <- function(x){
      return(evspline(B_f$os, sID = q, x = x)[,-1])
    }
    return(func)
  })
  
  cs = map(1:(Q-2), function(q){
    return(integral(function(x){return(x*basis[[q]](x))}, 
                    xmin = 0, xmax = 1))
  }) %>% unlist()
  slopeF = function(x){
    output = x - evspline(B_f$os, x = x)[,-1] %*% cs
    mag = 1/3 - sum(cs^2)
    return(output/sqrt(mag))
  }
  basis[[Q-1]] = slopeF
  
  ci = map(1:(Q-1), function(q){
    return(integral(basis[[q]], xmin = 0, xmax = 1))
  }) %>% unlist()
  intF = function(x){
    output = rep(1, length(x)) - evspline(B_f$os, x = x)[,-1] %*% ci[1:(Q-2)] - slopeF(x)*ci[Q-1]
    mag = 1 - sum(ci^2)
    return(output/sqrt(mag))
  }
  basis[[Q]] = intF
  return(list(B = basis, cInt = ci, cSlo = cs))
}

Splinet_d2 <- function(Q, cInt, cSlo){
  
  B_f = splinet(knots = seq(0, 10, length.out = Q+2), norm = T)
  second_deriv = deriva(deriva(B_f$os))
  
  derivs = map(1:(Q-2), function(q){
    func <- function(x){
      return(evspline(second_deriv, sID = q, x = x)[,-1])
    }
    return(func)
  })

  slopeD2 = function(x){
    output = - evspline(second_deriv, x = x)[,-1] %*% cSlo
    mag = 1/3 - sum(cSlo^2)
    return(output/sqrt(mag))
  }
  derivs[[Q-1]] = slopeD2
  
  intD2 = function(x){
    output = - evspline(second_deriv, x = x)[,-1] %*% cInt[1:(Q-2)] - slopeD2(x)*cInt[Q-1]
    mag = 1 - sum(cInt^2)
    return(output/sqrt(mag))
  }
  derivs[[Q]] = intD2
  return(derivs)
}
