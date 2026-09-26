% To test
% N=3;
% u=zeros(N+2,N+2);
% 
% u(3,2)=3.2;
% u(4,3)=4.3;
% 
% u
% 
% print_field(u);


function print_field(u)

[M,N]=size(u);

if (N<1000 && M<1000)
    for j = M:-1:1 % Va desde M y le va restando 1 hasta llegar a la pos. 1(se invierte la columna)
        fprintf('j=%2d  ',j);
        for i=1:N
            fprintf(' %+8.3e' , u(i,j));
        end
        fprintf('\n');
    end

end

end

% function test( u )
% print_field( u ) prints a field with the notation described in
% slides
%   u: field to print
% Example of use: print_field(up);
%     N=size(u,1)-2; % mesh size
%     for j=N+2:-1:1
%         fprintf('j=%2d  ',j);
%         for i=1:N+2
%             fprintf(' %+8.3e' ,u(i,j));
%         end
%         fprintf('\n');
%     end
% end
% 
