function computeUC_LSQ
    %%%%%%%%
    % Test least-squares algorithm for finding cell velocity from edge
    % values
    %%%%%%%%

    %%
    % Some suntans-like inputs
    
    %%%
    % Quad cell rotated at 45 degrees
    %%% 
    phys.u = [0.1,0.2,0.3,0.4]; % edge normal velocity
    grid.nfaces = 4;
    cf = 1/sqrt(2);
    grid.n1 = [cf,cf,-cf,-cf]; % normal of each face in x-direction
    grid.n2 = [cf,-cf,-cf,cf]; % normal of each face in y-direction

    %%%
    % Tri Cell
    %%%
%     phys.u = [0.1,0.2,0.3]; % edge normal velocity
%     grid.nfaces = 3;
%     grid.n1 = [0,0.5,0.5]; % normal of each face in x-direction
%     grid.n2 = [1,0.5,0.5]; % normal of each face in y-direction
    %%

    % Define the A.x = b components
    % A = [nface,2] matrix with normal components in each column
    % b = [nface,1] vector with edge normal velocity
    A = [grid.n1', grid.n2'];
    b = phys.u';

    %% Method 1) Use MATLAB backslash to solve least-squares problems
    uv_1 = A\b

    %% Method 2) Do the matrix transpose, then solve the linear system (again using backslash)
    % The least squares problem is solved via:
    %   A'x = b'
    % where,
    % A' = A^T * A
    % b' = A^T * b

    A_pr = A'*A;
    b_pr = A'*b;
    uv_2 = A_pr\b_pr

    %% Method 3) Same as method 2 but using the pseudo-c code matrix solver
    uv_3 = linsolve_lu(A_pr,b_pr,2)

    %% Method 4) Do the matrix multiplication by hand to be like c code
    % These operations can probably be done in place on the main vectors in the code??
    A_pr2 = zeros(2,2);
    b_pr2 = zeros(2,1);
    A_T = zeros(2,grid.nfaces);
    % A Matrix transpose
    for jj = 1:grid.nfaces
        for ii = 1:2
            A_T(ii,jj) = A(jj,ii);
        end
    end

    % Matrix multiplication
    for ii = 1:2
        for jj = 1:2
            ss=0;
            for kk = 1:grid.nfaces
                ss =  ss + A_T(ii,kk) * A(kk,jj);
            end
            A_pr2(ii,jj) = ss;

        end
    end

    for ii = 1:2
        for jj = 1:grid.nfaces
            b_pr2(ii) = b_pr2(ii) + A_T(ii,jj)*b(jj);
        end
    end

    uv_4 = linsolve_lu(A_pr2,b_pr2,2)
    
    %% Method 5) Do the matrix multiplication in place without initialising any intermediate arrays
    A_pr3 = zeros(2,2);
    b_pr3 = zeros(2,1);
    A_T = zeros(2,grid.nfaces);
    % A Matrix transpose
    for jj = 1:grid.nfaces
        for ii = 1:2
            A_T(ii,jj) = A(jj,ii);
        end
    end

    % Matrix multiplication
    for ii = 1:2
        for jj = 1:2
            ss=0;
            for kk = 1:grid.nfaces
                ss =  ss + A_T(ii,kk) * A(kk,jj);
            end
            A_pr3(ii,jj) = ss;

        end
    end

    for ii = 1:2
        for jj = 1:grid.nfaces
            b_pr3(ii) = b_pr3(ii) + A_T(ii,jj)*b(jj);
        end
    end

    uv_5 = linsolve_lu(A_pr3,b_pr3,2)

end

function b = linsolve_lu(A,b,N)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculates the solution,x, to A.x = b via LU decomposition
    % 
    % A is NxN square matrix
    % b is Nx1 vector
    % Returns Nx1 vector 'x'
    % Reference: Golub and Van Loan, Matrix Computations
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Algorithm to find lu decomp - See pa
    Atmp = A; % just overwrite A in c code
    for k = 1:N-1
        for i = k+1:N
            Atmp(i,k) = Atmp(i,k)/Atmp(k,k);
            for j = k+1:N
                Atmp(i,j) = Atmp(i,j)-Atmp(k,j)*Atmp(i,k);
            end
        end
    end

    L = zeros(N);
    U = zeros(N);
    for k = 1:N
        U(k,k)=Atmp(k,k);
        for j=k+1:N
            U(k,j)=Atmp(k,j);
        end
        L(k,k)=1.0;
        for j=1:k-1
            L(k,j)=Atmp(k,j);
        end
    end

    % Solve L.y=b via forward substitution (Alg 3.1.1);
    b(1) = b(1)/L(1,1);
    for i = 2:N
        sumi = 0;
        for j = 1:i-1
            sumi = sumi + L(i,j)*b(j);
        end
        b(i) = (b(i)-sumi) / L(i,i);
    end

    % Solve U.x=y via backward substitution (Alg 3.1.2)

    b(N) = b(N)/U(N,N);
    for i=N-1:-1:1
        sumi=0;
        for j=i+1:N
            sumi = sumi + U(i,j)*b(j);
        end
        b(i) = (b(i)-sumi)/U(i,i);
    end
    
end


