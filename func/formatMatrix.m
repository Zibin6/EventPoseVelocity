function str = formatMatrix(M)
    str = sprintf('[%.4f %.4f %.4f;\n\t             %.4f %.4f %.4f;\n\t             %.4f %.4f %.4f]', ...
    M(1,1), M(1,2), M(1,3), M(2,1), M(2,2), M(2,3), M(3,1), M(3,2), M(3,3));
end