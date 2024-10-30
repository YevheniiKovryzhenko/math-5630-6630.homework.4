hw04_worker = hw04();

hw_assert(abs(hw04_worker.p1([-1, 1;1, 1], 0) - 1) < 1e-8);
hw_assert(abs(hw04_worker.p1([-1, 1;1, 2], 0) - 1.5)<1e-8);
hw_assert(abs(hw04_worker.p1([-1, 1;1, 2;3, 3], 0) - 1.5)<1e-8);
hw_assert(abs(hw04_worker.p1([-1, 1;1, 2;3, 7], 0) - 1)<1e-8);
hw_assert(abs(hw04_worker.p1([-1, 1;1, 2;3, 7; 4, 2], 0) + 0.8)<1e-8);
hw_assert(abs(hw04_worker.p1([-1, 1;1, 2;3, 7; 4, 2; 5, 9], 0) + 6.25)<1e-8);


hw_assert(abs(hw04_worker.p2({[0, 0], [1,1,3,6]}, 0.5) - 1/8 ) < 1e-8);
hw_assert(abs(hw04_worker.p2({[0, 0], [1,1,4,12,24]}, 0.5) - 1/16 ) < 1e-8);
hw_assert((abs(hw04_worker.p2({[0, 0, 1], [1, 3, 6]}, [2, 3]) - [14, 39] ) < 1e-8) );



function hw_assert(X)
    if X; fprintf('\t PASS\n'); else; fprintf('\t FAIL\n'); end
end

