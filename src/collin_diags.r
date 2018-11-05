#------------
# Copyright (C) 2018 Brian J. Stucky.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
#------------


#------------
# Provides functions that implement various conditioning and collinearity
# diagnostics for linear models.  The methods are based on those developed in
# Belsley 1991.
#------------


#------
# Tests whether a matrix is singular.
#------
is_singular = function(matx) {
    return( class(try(solve(matx), silent=T)) == 'try-error' )
}


#------
# Scales a matrix so all columns have unit length under the 2-norm.
#------
scale_matrix = function(matx) {
  for (colnum in 1:ncol(matx)) {
    colvec = matx[,colnum]
    # If the norm of the vector is 0 (that is, the vector is a 0 vector), leave
    # it as is.
    if (sum(colvec) != 0) {
      matx[,colnum] = colvec / norm(colvec, type='2')
    }
  }
  
  return(matx)
}


#------
# Gets the condition number, condition indices, and variance decomposition
# proportions associated with a linear model formula and predictor variables.
#------
collin_diags = function(model_spec, data=environment(model_spec)) {
  d_mtx = scale_matrix(model.matrix(model_spec, data=data))
  
  xtx = t(d_mtx) %*% d_mtx
  if (is_singular(xtx)) {
    result = list(c_num=Inf, c_indices=NULL, VDPs=NULL)
    class(result) = c('CollinDiagsResults', class(result))
    return(result)
  }
  
  # Because we tested for the invertibility of X'X, X'X will not have 0 as an
  # eigenvalue, so all singular values of X will be nonzero.
  svdv = svd(d_mtx, nv=ncol(d_mtx))
  svals = svdv$d
  
  c_num = max(svals) / min(svals)
  
  # Calculate the condition indices.
  c_indices = max(svals) / svals
  
  # Calculate the variance decomposition proportions for each coefficient
  # in the regression model.
  var_decomp = matrix(nrow=0, ncol=ncol(d_mtx))
  for (varnum in 1:ncol(d_mtx)) {
    phis = (svdv$v[varnum,] ^ 2) / (svals ^ 2)
    var_decomp = rbind(var_decomp, phis / sum(phis))
  }
  colnames(var_decomp) = round(c_indices, 2)
  rownames(var_decomp) = colnames(d_mtx)

  result = list(
    c_num=c_num,
    c_indices=sort(c_indices, decreasing=T),
    VDPs=var_decomp
  )
  class(result) = c('CollinDiagsResults', class(result))
  
  return(result)
}

#------
# Defines a custom print function for collinearity diagnostics result objects.
#------
print.CollinDiagsResults = function(result) {
  cat('Condition number of design matrix: ', result$c_num, '\n\n', sep='')
  cat('Condition indices:\n')
  cat(result$c_indices, sep=', ')
  cat('\n\nVariance Decomposition Proportions:\n')
  print(round(result$VDPs, 2))
}

