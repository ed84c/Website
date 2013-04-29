function p2Move = NandCver5Player2Strategy(board)
board(board==1)=1;
board(board==-1)=-1;
availableMoves = zeros(1,9);
tempP2Move = 0;

%Find Available Moves
availableMoves = findAvailable(board)

%TryToWin
tempP2Move = winnerMove(board,-1)

if tempP2Move == 0
	printf('Attempt block')
	tempP2Move = winnerMove(board,1)
end

%Branch Case 1
if tempP2Move == 0 
	printf('Branch Case 1 \n')
d1 = diag(board);
d2 = diag(rot90(board));
	if isequal(d1,[1 -1 1]')
	availableMoves(availableMoves==3) = [];
	availableMoves(availableMoves==7) = [];
	end
	if isequal(d2,[1 -1 1]')
	availableMoves(availableMoves==1) = [];
	availableMoves(availableMoves==9) = [];
	end
end
availableMoves
%Branch Case 2
if tempP2Move ==0
	printf('Branch Case 2 \n')
	tempP2Move = branchBlock(board,availableMoves)
end

%Go centre square
if tempP2Move == 0 && board(2,2)==0
	printf('Go centre')
	tempP2Move = 5;
end

%Go anywhere
if tempP2Move == 0
	printf('Go anyhwhere')
	tempP2Move = availableMoves(1);
end

p2Move = tempP2Move;

end

function availableMoves = findAvailable(board);
linBoard = zeros(1,9);
linBoard = reshape(board,1,9);
moveVector = [1,2,3,4,5,6,7,8,9]
linBoard = ~logical(linBoard);
availableMoves = moveVector.*linBoard;
availableMoves(availableMoves == 0)=[];
end


function winningMove = winnerMove(board,defensive)
tempBoard = zeros(3,3);
tempBoard = board;
i=0;
line1 = zeros(3,1);
line2 = zeros(3,1);
line3 = zeros(3,1);
sum1 =0;
sum2=0;
sum3=0;
emptyBoard = zeros(3,3);
linEmptyBoard = zeros(1,9);
currMove = 0;

while i<6 && sum1*defensive<2 && sum2*defensive<2 && sum3*defensive<2
	line1= tempBoard(:,1);
	sum1 = sum(line1);
	
	line2 = tempBoard(:,2);
	sum2 = sum(line2);
	
	line3 = diag(tempBoard);
	sum3 = sum(line3);
	
	i=i+1;
	tempBoard = rot90(tempBoard);
end
sums = [sum1 sum2 sum3]
if 	sum1*defensive ==2
	printf('Outer Line\n')
	line1(line1==0)=3;
	emptyBoard(:,1) = line1;
	emptyBoard = rot90(emptyBoard,(i-1)*(-1))
	linEmptyBoard = reshape(emptyBoard,1,9)
	currMove = find(linEmptyBoard==3)
end

if 	sum2*defensive ==2
	printf('Inner Line\n')
	line2(line2==0)=3;
	emptyBoard(:,2) = line2;
	emptyBoard = rot90(emptyBoard,(i-1)*(-1))
	linEmptyBoard = reshape(emptyBoard,1,9)
	currMove = find(linEmptyBoard==3)
end

if 	 sum3*defensive ==2
	printf('Diagonal')
	line3(line3==0)=3;
	emptyBoard = blkdiag(line3(1),line3(2),line3(3))
	emptyBoard = rot90(emptyBoard,(i-1)*(-1));
	linEmptyBoard = reshape(emptyBoard,1,9);
	currMove = find(linEmptyBoard==3)
end
	winningMove = currMove
end

function branchBlockMove = branchBlock(board, availableMoves)
tempBoard = board;
i=0;
bBranch = 0;
currMove = 0;

while i<6 && bBranch ==0
	if sum(tempBoard(:,1)) == 1 && sum(tempBoard(1,:))==1 && tempBoard(1,1) ==0
		bBranch=1
	end
	i=i+1;
	tempBoard = rot90(tempBoard);
end

if bBranch == 1
	switch(i)
	case 1
		currMove = 1;
	case 2
		currMove =7;
	case 3
		currMove = 9;
	case 4
		currMove = 3;
	otherwise
		currMove = 0;
	end
end
	if find(availableMoves==currMove) >0
		branchBlockMove = currMove;
	else
		branchBlockMove = 0
	end
end