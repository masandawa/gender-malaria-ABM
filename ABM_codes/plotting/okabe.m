function c = okabe(name)
%OKABE  Colour-blind-safe Okabe-Ito colours.
switch name
    case 'female', c = [213 94 0]/255;     % vermillion
    case 'male',   c = [0 114 178]/255;    % blue
    case 'preg',   c = [204 121 167]/255;  % reddish purple
    case 'green',  c = [0 158 115]/255;
    case 'orange', c = [230 159 0]/255;
    case 'sky',    c = [86 180 233]/255;
    case 'grey',   c = [0.55 0.55 0.55];
    case 'band',   c = [0.88 0.94 0.90];
    otherwise,     c = [0 0 0];
end
end
