local ok, mtoc = pcall(require, "mtoc")
if ok then
    mtoc.update_config({
        toc_list = {
            markers = "1.",
            indent_size = 3,
        },
    })
end
