function NandC
warning('off','all');
N = input('Number of games to play');
wins = zeros(2,1);
winsPerc = zeros(2,1);

for i=1:N
	winner = playgame;
	switch (winner)
		case(0)
		case(1)
			wins(1) = wins(1) + 1;
			break
		case(2)
			wins(2) = wins(2) + 1;
	end
end

winsPerc = wins/N * 100;

printf('%i games were played with: \n player 1 winning %i (%-2.2f %%) \n player 2 winning %i (%-2.2f %%) \n Draws: %i (%-2.2f%%) \n'...
,N,wins(1),winsPerc(1),wins(2),winsPerc(2),N-wins(1)-wins(2),100-winsPerc(1)-winsPerc(2));

			
end

			
function [winner] = playgame()
M = zeros(3,3);
Mlin = zeros(9,1);
win = 0;
i=0;

while (i<9&&win==0)

	if (mod(i,2)==0)
		Mlin = reshape(M,9,1);
		player1Move = player1(M);
		if Mlin(player1Move)~=0;
			printf('Error')
			break;
		else
			Mlin(player1Move) = 1;
			currPlayer =1;
		end
	else
		Mlin = reshape(M,9,1);		
		player2Move = player2(M)
		
		if Mlin(player2Move)~=0
			printf('Error')
			break;
		else
			Mlin(player2Move) = -1;
			currPlayer=2;
		end
	end
	
	M = vec2mat(Mlin,3);
	M=M';
	%pause(0.5);
	i=i+1;
	win = winchecker(M);
	%Mhistory(i) = M;
	outputBoard(M)
end

if (win==true)
	%printf('The winner is %i \n', currPlayer);
	winner = currPlayer;
else
	%printf('The game was a draw \n');
	winner = 0;
end


end


function [move_p1] = player1(board)
move_p1_index=0;
available_movesp1 = 0;
available_movesp1 = zeros(9);
available_movesp1 = available_moves(board);
n_available_moves = size(available_movesp1,1);

move_p1_index = randi([1 n_available_moves],1,1);
move_p1=available_movesp1(move_p1_index);

end

function [move_p2] = player2(board)
available_movesp2 = 0;
move_p2_index=0;
available_movesp2 = zeros(9);
available_movesp2 = available_moves(board);
n_available_moves = size(available_movesp2,1);

a =0;
tempmove_p2 = 0;

tempmove_p2 = NandC_Player2Strategy(board, -1);

if tempmove_p2 ==0
	tempmove_p2 = NandC_Player2Strategy(board, 1);
end 

if (tempmove_p2==0&&board(2,2)==0)
	tempmove_p2 = 5;
end

freeCornerIndices = zeros(1,4);
freeCornerIndicies = find(available_movesp2==1 | available_movesp2==3 | available_movesp2==7 | available_movesp2==9)
if (tempmove_p2==0 & size(freeCornerIndicies)>0)
	tempmove_p2 = available_movesp2(freeCornerIndicies (1));
end

boardTemp = board;
branchFound = 0;
cornerPositions = [1 3 9 7]
i=0;
while (i<5 && branchFound == 0)
	cornerValues = [boardTemp(1), boardTemp(2), boardTemp(3), boardTemp(4), boardTemp(5)]
	boardTemp = rot90(boardTemp);
	if (sum(cornerValues) == 2 && find(available_movesp2 == cornerPositions(i)) ~=[]);
		tempmove_p2 = cornerPoistions(i);
		branchFound = 1;
	end
end


if diag(board) == [-1 1 -1]
	available_movesp2(avaialable_movesp2==3) = [];
	available_movesp2(available_movesp2==7)=[];
	n_available_moves = size(available_movesp2,1);
end

if diag(rot90(board)) == [-1 1 -1]
	avaiable_movesp2(available_movesp2==1) = [];
	available_movesp2(available_movesp2==9) = [];
	n_available_moves = size(available_movesp2,1);	
end 

if tempmove_p2==0
	move_p2_index = randi([1 n_available_moves],1,1);
	tempmove_p2=available_movesp2(move_p2_index);
end 

	move_p2 = tempmove_p2;
end

function [available_move] = available_moves(board)
linboard= reshape(board,9,1);
B = [1;2;3;4;5;6;7;8;9];
available_moveTemp = ~(logical(linboard)).*B;
available_moveTemp(available_moveTemp==0)=[];
available_move = available_moveTemp;
end

function [wincheck]  = winchecker(board)

Atemp = int8(board);
i=0;
sumColOuter=0;
sumColInner=0;
sumColDiag=0;

while (i<5 & abs(sumColOuter)<3 & abs(sumColInner)<3 & abs(sumColDiag)<3)
	sumColOuter = sum(Atemp(1,:));
	sumColInner = sum(Atemp(2,:));
	sumColDiag =  trace(Atemp);
	Atemp = rot90(Atemp);
	i=i+1;
end

if (i<5)
	wincheck = 1;
else
	wincheck = 0;
end

end  

function outputBoard(board)
drawnBoard = zeros(3,3);
drawnBoard = board;
drawnBoard(drawnBoard==1)='x';
drawnBoard(drawnBoard==-1)='O';
drawnBoard(drawnBoard==0)=' ';

printf('%c|%c|%c',drawnBoard(1,1), drawnBoard(1,2), drawnBoard(1,3));
printf('\n-----\n');
printf('%c|%c|%c',drawnBoard(2,1), drawnBoard(2,2), drawnBoard(2,3));
printf('\n-----\n');
printf('%c|%c|%c',drawnBoard(3,1), drawnBoard(3,2), drawnBoard(3,3));
printf('\n \n \n');


end

function soundP1=soundP1;
x=[0:1/20000:0.1];
soundP1 = sin(x*5000);
end 

function soundP2=soundP2;
x=[0:1/20000:0.1];
soundP2 = sin(x*5000);
end 