% ROUND TO THE NEAREST EVEN INTEGER

% COPYRIGHT: 
%
%   MIT License
%   © 2023 Dario Chiantello <pipesoundlibrary@gmail.com>
%

function s = roundEven(s)

mod2 = mod(s,2);
s = s - mod2;
