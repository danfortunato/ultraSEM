function T = alphabet( letter )
% ALPHABET   Define a letter-shaped domain.
%
% T = ALPHABET(LETTER) defines an ULTRASEMDOMAIN object that represents the
% letter in LETTER. Currently, we only support the 26 letters in the Latin
% alphabet.

switch letter

    case '-'
        dom = [-1 0 0 1 ;
            0 1 0 1 ];
        idx = {[1 2]};

    case '!'
    dom = [ 0 1 0 1 ;
            0 1 1 2;
            0 1 2 3
            0 1 -2 -1];
        idx = {[1 2 ; 3 NaN ; 4 NaN], ...
            [1 2 ; 3 NaN], [1 2]};

    case 'A'
        dom = [0 1 0 1 ;
            0 1 1 2 ;
            0 1 2 3 ;
            2 3 0 1 ;
            2 3 1 2 ;
            2 3 2 3 ;
            1 2 1 2 ;
            1 2 2 3];
        idx = {[1 2 ; 3 NaN ; 4 5 ; 6 NaN ; 7 8], ...
            [1 2 ; 3 4 ; 5 NaN ], [1 2 ; 3 NaN], [1 2]};
           
    case 'E'
        dom = [0 1 0 1 ;
            0 1 1 2 ;
            0 1 2 3 ;
            0 1 3 4 ; 
            0 1 4 5 ; 
            1 2 0 1  ;
            2 3 0 1 ; 
            1 2 2 3 
            1 2 4 5 
            2 3 4 5 
            ];
        idx = {[1 2 ; 3 8 ; 4 5 ; 6 7 ; 9 10], ...
            [1 4 ; 2 NaN ; 3 5 ], [1 2 ; 3 NaN], [1 2]};           

    case 'H'
        dom = [0 1 0 1 ;
            0 1 1 2 ;
            0 1 2 3 ;
            2 3 0 1 ;
            2 3 1 2 ;
            2 3 2 3 ;
            1 2 1 2 ];
        idx = {[1 2 ; 3 NaN ; 4 5 ; 6 NaN ; 7 NaN], ...
            [1 2 ; 3 4 ; 5 NaN ], [1 2 ; 3 NaN], [1 2]};

    case 'l'
        dom = [-1 0 -1 0 ;
            0 1 -1 0 ;
            -1 0  0 1];
        idx = {[1 2 ; 3 NaN], ...
            [1 2]};

    case 'L'
        dom = [-1 0 -1 0 ;
            0 1 -1 0 ;
            -1 0 0 1 ;
            -1 0 1 2 ];
        idx = {[1 2 ; 3 4], ...
            [1 2]};

    case 'M'
        dom = [0 1 0 1 ;
            0 1 1 2 ;
            0 1 2 3 ;
            0 1 3 4 ;
            0 1 4 5 ;
            2 3 0 1 ;
            2 3 1 2 ;
            2 3 2 3 ;
            2 3 3 4 ;
            2 3 4 5 ;            
            1 2 2 3 
            1 2 3 4 ];
        idx = {};
        idx = {[1 2 ; 3 NaN ; 4 5 ; 6 7 ; 8 NaN ; 9 10 ; 11 12], ...
            [1 2 ; 3 7 ; 4 5 ; 6 NaN ], [1 2 ; 3 4], [1 2]};     

    case 'o'
        dom = [-1 0 -1 0
            0 1 -1 0
            1 2 -1 0
            1 2 0 1
            1 2 1 2
            0 1 1 2
            -1 0 1 2
            -1 0 0 1];
        idx = {[1 2 ; 3 4 ; 5 6 ; 7 8], ...
            [1 2 ;  3 4], ...
            [1 2]};

    case 'O'

        dom = [-2 -1 -1 0
            -2 -1 -2 -1
            -1 0 -2 -1
            0 1 -2 -1
            1 2 -2 -1
            1 2 -1 0
            1 2 0 1
            1 2 1 2
            0 1 1 2
            -1 0 1 2
            -2 -1 1 2
            -2 -1 0 1];
        idx = {[1 2 ; 3 NaN ; 4 NaN ; 5 6 ; 7 8 ; 9 NaN ; 10 NaN ; 11 12], ...
            [1 2 ; 3 4 ; 5 6; 7 8], ...
            [1 2 ; 3 4], ...
            [1 2]};

    case 'p'
        dom = [-1 0 -1 0
            0 1 -1 0
            1 2 -1 0
            1 2 0 1
            1 2 1 2
            0 1 1 2
            -1 0 1 2
            -1 0 0 1
            -1 0 -2 -1];
        idx = {[1 2 ; 3 4 ; 5 6 ; 7 8 ; 9 NaN], ...
            [1 2 ;  3 4 ; 5 NaN], ...
            [1 2 ; 3 NaN], ...
            [1 2]};

    case 'P'
        dom = [-1 0 -1 0
            0 1 -1 0
            1 2 -1 0
            1 2 0 1
            1 2 1 2
            0 1 1 2
            -1 0 1 2
            -1 0 0 1
            -1 0 -2 -1
            -1 0 -3 -2];
        idx = {[1 2 ; 3 4 ; 5 6 ; 7 8 ; 9 10], ...
            [1 2 ;  3 4 ; 5 NaN], ...
            [1 2 ; 3 NaN], ...
            [1 2]};

    case {'s'}
        dom = [0 1 0 1 ;
            1 2 0 1 ;
            1 2 1 2 ;
            1 2 2 3 ;
            0 1 2 3 ;
            0 1 3 4 ;
            0 1 4 5 ;
            1 2 4 5 ];
        idx = {[1 2 ; 3 4 ; 5 6 ; 7 8], ...
            [1 2 ; 3 4], ...
            [1 2]};
        
    case {'S'}
        dom = [0 1 0 1 ;
            1 2 0 1 ;
            2 3 0 1 ;
            2 3 1 2 ;
            2 3 2 3 ;
            1 2 2 3 ;
            0 1 2 3 ;
            0 1 3 4 ;
            0 1 4 5 ;
            1 2 4 5 
            2 3 4 5 ];
        idx = {[1 2 ; 3 4 ; 5 6 ; 7 8 ; 9 10 ; 11 NaN], ...
            [1 2 ; 3 4 ; 5 6], ...
            [1 2 ; 3 NaN], [ 1 2]};        

    case 't'
        dom = [0 1 1 2 ;
            1 2 1 2 ;
            2 3 1 2 ;
            1 2 0 1];
        idx = {[1 2 ; 3 NaN ; 4 NaN], ...
            [1 2 ; 3 3], ...
            [1 2]};

    case 'T'
        dom = [0 1 1 2 ;
            1 2 1 2 ;
            2 3 1 2 ;
            1 2 0 1 ;
            1 2 -1 0];
        idx = {[1 2 ; 3 NaN ; 4 5], ...
            [1 2 ; 3 NaN], ...
            [1 2]};

    case 'u'
        dom = [0 1 0 1;
            1 2 0 1 ;
            2 3 0 1 ;
            2 3 1 2 ;
            2 3 2 3 ;
            0 1 1 2 ;
            0 1 2 3 ];
        idx = {[1 2 ; 3 4 ; 5 NaN ; 6 7], ...
            [1 4 ; 2 3], ...
            [1 2]};

    case 'U'
        dom = [0 1 0 1;
            1 2 0 1 ;
            2 3 0 1 ;
            2 3 1 2 ;
            2 3 2 3 ;
            2 3 3 4 ;
            2 3 4 5 ;            
            0 1 1 2 ;
            0 1 2 3 ;
            0 1 3 4
            0 1 4 5];
        idx = {[1 2 ; 3 4 ; 5 6 ; 7 NaN ; 8 9 ; 10 11], ...
            [1 2 ; 3 4 ; 5 6], ...
            [1 2 ; 3 NaN], [1 2]};

    otherwise
        error('Unknown letter');

end

dom = ultraSEMRect(dom);
T = ultraSEMDomain(dom, idx);

end
