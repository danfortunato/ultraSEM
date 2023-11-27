function D = L()

D1 = ultraSEM.rectangle([-1 0 -1 0]);
D2 = ultraSEM.rectangle([0 1 -1 0]);
D3 = ultraSEM.rectangle([-1 0 0 1]);
D = (D1&D2) & D3;

end