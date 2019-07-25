%%%%%% GETTING STARTED

% The first thing you should do is place this file in a folder that you
% will use to save your work.

% Now run Matlab and use the "Current Directory" pull-down menu to point 
% Matlab to the folder where you have saved this file.

% You will need to work with two pieces of Matlab, the Editor and the
% Command Window.

% The Command Window should open automatically when you run Matlab.

% To open the Editor, type "edit" and hit return at the command prompt.

% You use the Editor to create, open, save, etc, Matlab program files.
% These files have the extension ".m"

% If you have a program file called "my_program.m", you run it by typing
% "my_program" and then hitting enter in the Command Window. Of course,
% Matlab has to be pointed to the right folder on your computer (ie, the
% folder where "my_program.m" is saved).

% To see a list of all the files in the directory, type "dir" at the
% Command Window.

% To see a list of all the Matlab objects in the current directory, type
% "who" or "whos". "who" gives a list of all the named variables in the
% current directory. "whos" gives a list with details, eg, the size of
% various matrices, how much memory they take to store, and what kind of
% object they are (eg, a symbol).

% In order for there to be any Matlab objects, we need to create some.

%%%%% MATRICES AND VECTORS

% We can create a 2-by-2 matrix in the following way

A = [1 2;3 4]

% When you run this program, you will see the matrix A appear in the
% command window. Bigger matrices can be made in an analogous way. 

% The semi-colon starts a new row. Here is another matrix

B = [5 6;7 8];

% Putting a semi-colon at the end of a line stops Matlab from displaying
% the answer in the command window. 

% Adding and subtracting are easy

C = A + B
D = A - B

% Of course, matrices must be conformable or Matlab will give an error.

% A prime (') indicates transposition

Ctranspose  = C'

% This also turns column vectors into row vectors

column = [1;2;3;4]

row = column'

% Multiplying is also easy

E = A*B;

% Whenever you want to see something in the command window, delete the
% semi-colon at the end of the line, save this file again, and run it.

% Alternatively, you could just type "E" (or "A" or whatever the object 
% you're interested in is) at the command line one the
% program is finished and Matlab will call up the object from its memory.

% Taking a matrix inverse works so long as the matrix is not singular.
% Compare

invC = inv(C)

% with 

invD = inv(D)

% (this last command gives an error message since D is singular). To see
% check that D is singular, note

detD = det(D)

% is zero

% Matrices can also be created with a number of specialized commands. For
% example:

O = ones(3,3);

% creates a 3-by-3 matrix with all elements equal to one, 

Z = zeros(3,3);

% creates a 3-by-3 matrix of zeros, and

I = eye(3,3);

% creates a 3-by-3 identity matrix 

% We can also stack matrices together. Say,

bigmatrix = [O Z;Z O];

% We need to make sure that all our row and column dimensions add up. It is
% sometimes useful to check things like

[m,n] = size(bigmatrix);

% which gives the row (=m) and column (=n) dimensions.

%%%%%% OPERATIONS

% We can do commands like max,min,sum and apply them to matrices. These are
% done column by column. 

maxA = max(A);
minA = min(A);
sumA = sum(A);

% Some other helpful operations are cumulative sums, eg

time = cumsum(ones(10,1));

% creates a vector that runs from 1 to 10.

%%%%% LOOPS

% The basic "for" loop has the structure

N  = 10;         %% some number
cs = zeros(N,1); %% a vector to store some answers

for i = 1:N
    cs(i) = sum(ones(i,1));
end

% This produces a vector that is the same as cumsum(ones(10,1))

% The basic "while" loop has the structure

tolerance = 10^-6;  %% an arbitrary tolerance level
error     = 100;    %% some number bigger than tolerance, to get us going
n         = 1;      %% an initial condition
xs        = [];     %% an empty matrix -- we will fill this in as we go

while error > tolerance
    xs_new = 0.5^n;   %% calculate some value, say 0.5^n -- a convergent series
	xs = [xs;xs_new]; %% keep adding new elements to xs
    
    error = abs(xs_new-0); %% check to see if we have convergence
    n     = n+1;           %% if not, keep iterating!
    
end

%%%%% PLOTS

% We can see what this sequence looks like by

plot(xs)
title('0.5^n')
xlabel('n')













