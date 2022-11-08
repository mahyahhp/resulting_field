%% 
% === inputs ================================================================
% "threshold", "c", "d":    A variable 
% "F":                      A scalar field - the size of this matrix is N x M
% "B":                      A scalar field - the size of this matrix is N x M
% "R":                      A scalar field - the size of this matrix is N x M
% ===========================================================================
% 
% === output ================================================================
% "R":                      A scalar field - the size of this matrix is N x M
% ===========================================================================

function [R] = field_calculator(threshold, ... 	% in
								c, ... 			% in
								d, ... 			% in
								F, ... 			% in
								B, ... 			% in
                                R)     			% in-out  

	%% --- Variables needed for the calculations ----------------------------------------------------------
	error_tol=1e-12; % the tolerance value needed in equality conditions
	% ALTERNATIVE METHOD: This variable (error_tol) can be an input to the function 
	%                     because its value depends on the considered decimal precision of the solution 
	%                     (here a double precision is assumed) and might be the same throughout the software. 
	
	m_plus_11 = -2; m_plus_12 = -1; m_plus_13 = -2; 
	m_plus_21 = -1; m_plus_22 = 12; m_plus_23 = -1; 
	m_plus_31 = -2; m_plus_32 = -1; m_plus_33 = -2;
	
	m_minus_11 =  0;  m_minus_12 = -1;  m_minus_13 =  0; 
	m_minus_21 = -1;  m_minus_22 =  4;  m_minus_23 = -1; 
	m_minus_31 =  0;  m_minus_32 = -1;  m_minus_33 =  0;
	
	%% --- Resulting field calculations -------------------------------------------------------------------
	[N,M]=size(F);
	% ALTERNATIVE METHOD: If N and M are known, then these two can come as inputs to the function
	%                     without the need of redoing anything like
    %                     [N,M]=size(F). So, a faster code!

	parfor i=1:N 
		for j=1:M

			% PERFORMACE BOTTLENECK: If the matrix size is not large enough,
			%                        then the parallel overhead may dominate calculations, so the
			%                        expected scalability of the code would not be achieved, or even the
			%                        serial code might perform better than the parallel code.
	
			if (abs(B(i,j)-c)<=error_tol)
				R(i,j)=B(i,j);
				continue;
			elseif (abs(B(i,j)-d)<=error_tol)
				continue;
			end
			
			if((F(i,j)-B(i,j))>=threshold)
				if(i~=1 && i~=N && j~=1 && j~=M)
					R(i,j)=F(i-1, j-1) * m_plus_11 + F(i-1, j) * m_plus_12 + F(i-1, j+1) * m_plus_13 + ...
					       F(i  , j-1) * m_plus_21 + F(i  , j) * m_plus_22 + F(i  , j+1) * m_plus_23 + ...
					       F(i+1, j-1) * m_plus_31 + F(i+1, j) * m_plus_32 + F(i+1, j+1) * m_plus_33;
				else                
					R(i,j)=F(i,j);
					% NOTE: For boundary cells, R(i,j)=F(i,j) is assumed, because there is no
					%       a mask/stencil matrix designed for these boundary cells. In
                    %       terms of performance, it is also preferred to
                    %       separate the loop of boundary cells to avoid
                    %       extra "if statements" in this parallel loop
                    %       that might slowdown the computations.
                    %       However, as requested a short code in the test definition,
                    %       the boundary cells are handled here, although it is not efficient. 
				end
            else %(F(i,j)-B(i,j))<threshold
				if(i~=1 && i~=N && j~=1 && j~=M)
					R(i,j)=F(i-1, j-1) * m_minus_11 + F(i-1, j) * m_minus_12 + F(i-1, j+1) * m_minus_13 + ...
					       F(i  , j-1) * m_minus_21 + F(i  , j) * m_minus_22 + F(i  , j+1) * m_minus_23 + ...
					       F(i+1, j-1) * m_minus_31 + F(i+1, j) * m_minus_32 + F(i+1, j+1) * m_minus_33;
				else
					R(i,j)=F(i, j);
				end
			end
			% PERFORMACE BOTTLENECK: Some load imbalance! A worker,
			%                        may do nothing, depending on a satisfied condition
			%                        (for example if "abs(B(i,j)-d)<=error_tol"), 
			%                        but another worker may do a lot 
			%                        (for example if "(F(i,j)-B(i,j))>=threshold" and "i~=1 && i~=N && j~=1 && j~=M)" 
		end
    end

end



