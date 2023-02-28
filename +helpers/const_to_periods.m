function T= const_to_periods(constituents)
            T = [];
            for i = 1:numel(constituents)
                name = constituents{i};
                switch name
                    case 'M2'
                        t = 12.4206012;
                    case 'S2'
                        t = 12;
                    case 'N2'
                        t = 12.65834751;
                    case 'K1'
                        t = 23.93447213;
                    case 'M4'
                        t = 6.210300601;
                    case 'O1'
                        t = 25.81933871;
                    case 'M6'
                        t = 4.140200401;
                    case 'MK3'
                        t = 8.177140247;
                    case 'S4'
                        t = 6;
                    case 'MN4'
                        t = 6.269173724;
                    otherwise
                        error('Unknown tidal constituent')
                end
                T(i) = t;
            end