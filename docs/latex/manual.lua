-- pandoc filter for the GZero manual:
-- * the first "# Title" and the line after it become the title page
-- * ```latex blocks become display math
-- * table columns get widths from their content, so every table fills the text width
-- * long rows in tables wrap instead of running off the page

local title_done = false

-- inline code: allow line breaks after _ / . so long file names wrap inside table cells
local special = { ["\\"] = "\\textbackslash{}", ["{"] = "\\{", ["}"] = "\\}", ["$"] = "\\$",
  ["&"] = "\\&", ["#"] = "\\#", ["^"] = "\\textasciicircum{}", ["~"] = "\\textasciitilde{}",
  ["%"] = "\\%", ["<"] = "\\textless{}", [">"] = "\\textgreater{}", ["|"] = "\\textbar{}",
  ["_"] = "\\_\\allowbreak{}", ["/"] = "/\\allowbreak{}", ["."] = ".\\allowbreak{}" }

local function code_to_latex(el)
  local out = el.text:gsub("[\\{}$&#^~%%<>|_/.]", special)
  return pandoc.RawInline("latex", "\\texttt{" .. out .. "}")
end


function Pandoc(doc)
  local blocks = doc.blocks
  if blocks[1] and blocks[1].t == "Header" and blocks[1].level == 1 then
    doc.meta.title = pandoc.MetaInlines(blocks[1].content)
    table.remove(blocks, 1)
    if blocks[1] and blocks[1].t == "Para" then
      doc.meta.subtitle = pandoc.MetaInlines(blocks[1].content)
      table.remove(blocks, 1)
    end
  end
  -- after the tables are measured: quantities to math, inline code to breakable LaTeX
  doc.blocks = doc.blocks:walk(MATH_WALK)
  doc.blocks = doc.blocks:walk({ Code = code_to_latex })
  -- "## 1. ..." is a chapter: one level up (section), "###" a subsection
  doc.blocks = doc.blocks:walk({
    Header = function(h)
      h.level = math.max(1, h.level - 1)
      return h
    end
  })
  return doc
end

function CodeBlock(el)
  if el.classes[1] == "latex" then
    return pandoc.RawBlock("latex", "\\[\n" .. el.text .. "\n\\]")
  end
end

local function cell_len(cell)
  return #pandoc.utils.stringify(cell.contents)
end

function Table(tbl)
  local ncol = #tbl.colspecs
  local maxlen, maxtok = {}, {}
  for i = 1, ncol do maxlen[i] = 4; maxtok[i] = 3 end
  local function scan(rows)
    for _, row in ipairs(rows) do
      for i, cell in ipairs(row.cells) do
        if i <= ncol then
          local text = pandoc.utils.stringify(cell.contents)
          maxlen[i] = math.max(maxlen[i], #text)
          -- longest piece that cannot be broken (code breaks after _ / .)
          for tok in text:gmatch("[^%s_/%.]+[_/%.]?") do maxtok[i] = math.max(maxtok[i], #tok) end
        end
      end
    end
  end
  scan(tbl.head.rows)
  for _, body in ipairs(tbl.bodies) do scan(body.body) end
  -- square root keeps short columns readable and long ones from taking everything;
  -- every column is at least as wide as its longest unbreakable word (~0.0135 of the width per character)
  local w, total = {}, 0
  for i = 1, ncol do w[i] = math.sqrt(math.min(maxlen[i], 400)); total = total + w[i] end
  local frac, floor_sum, free = {}, 0, 0
  for i = 1, ncol do
    local f = 0.97 * w[i] / total
    local fmin = math.min(0.6, 0.0135 * maxtok[i] + 0.02)
    if f < fmin then frac[i] = fmin; floor_sum = floor_sum + fmin else frac[i] = nil; free = free + w[i] end
  end
  for i = 1, ncol do
    if frac[i] == nil then frac[i] = (0.97 - floor_sum) * w[i] / free end
    tbl.colspecs[i] = {pandoc.AlignLeft, frac[i]}
  end
  return tbl
end

-- plain-text quantities in running text (eps0, Gamma_L, G0, V_t+, x^2 ...) become math in the PDF.
-- The markdown stays plain so it reads well on GitHub. Italic text (setting names as shown in the
-- program) and code are left alone.
local SUBS = {
  {"%f[%w]Gamma_L,R%f[^%w_]", "\\Gamma_{L,R}"},
  {"%f[%w]Gamma_tot%f[^%w_]", "\\Gamma_\\mathrm{tot}"},
  {"%f[%w]Gamma_L%f[^%w_]", "\\Gamma_L"}, {"%f[%w]Gamma_R%f[^%w_]", "\\Gamma_R"},
  {"%f[%w]eps0%^2%f[^%w_]", "\\varepsilon_0^2"}, {"%f[%w]Gamma%^2%f[^%w_]", "\\Gamma^2"},
  {"%f[%w]eps0%f[^%w_]", "\\varepsilon_0"}, {"%f[%w]Gamma%f[^%w_]", "\\Gamma"},
  {"%f[%w]V_trans%+", "V_\\mathrm{trans}^{+}"}, {"%f[%w]V_trans%-", "V_\\mathrm{trans}^{-}"},
  {"%f[%w]V_trans%f[^%w_]", "V_\\mathrm{trans}"},
  {"%f[%w]V_t%+", "V_t^{+}"}, {"%f[%w]V_t%-", "V_t^{-}"}, {"%f[%w]V_t%f[^%w_]", "V_t"},
  {"%f[%w]sigma_log%f[^%w_]", "\\sigma_\\mathrm{log}"},
  {"%f[%w]G0%f[^%w_]", "G_0"}, {"%f[%w]J0%f[^%w_]", "J_0"}, {"%f[%w]n_j%f[^%w_]", "n_j"},
  {"%f[%w]10%^mu%f[^%w_]", "10^{\\mu}"}, {"%f[%w]f%^alpha%f[^%w_]", "f^{\\alpha}"},
  {"%f[%w]beta%f[^%w_]", "\\beta"}, {"%f[%w]sigma%f[^%w_]", "\\sigma"}, {"%f[%w]kappa%f[^%w_]", "\\kappa"},
  {"%f[%w]alpha%f[^%w_]", "\\alpha"}, {"%f[%w]mu%f[^%w_]", "\\mu"}, {"%f[%w]tau%f[^%w_]", "\\tau"},
  {"%f[%w]pi%f[^%w_]", "\\pi"},
  {"%f[%w](%d+)%^(%-?%d+)", "%1^{%2}"},
  {"%f[%w]([A-Za-z])%^(%w+)", "%1^{%2}"},
}
local TEXT_SUBS = { {"%f[%w]cm%^2%f[^%w_]", "cm$^2$"}, {"%f[%w]cm%^%-2%f[^%w_]", "cm$^{-2}$"} }

local function tex_escape(s)
  return (s:gsub("[\\{}$&#^_~%%]", { ["\\"] = "\\textbackslash{}", ["{"] = "\\{", ["}"] = "\\}",
    ["$"] = "\\$", ["&"] = "\\&", ["#"] = "\\#", ["^"] = "\\textasciicircum{}", ["_"] = "\\_",
    ["~"] = "\\textasciitilde{}", ["%"] = "\\%" }))
end

local function math_str(el)
  local text, store = el.text, {}
  local function keep(latex) store[#store + 1] = latex; return "\1" .. #store .. "\2" end
  for _, r in ipairs(TEXT_SUBS) do
    text = text:gsub(r[1], function(...) return keep((r[2]:gsub("%%(%d)", {...}))) end)
  end
  for _, r in ipairs(SUBS) do
    text = text:gsub(r[1], function(...)
      local caps = {...}
      local out = r[2]:gsub("%%(%d)", function(k) return caps[tonumber(k)] end)
      return keep("$" .. out .. "$")
    end)
  end
  if #store == 0 then return nil end
  local latex = tex_escape(text):gsub("\1(%d+)\2", function(k) return store[tonumber(k)] end)
  return pandoc.RawInline("latex", latex)
end

MATH_WALK = {
  traverse = "topdown",
  Emph = function(e) return e, false end,
  Code = function(e) return e, false end,
  Link = function(e) return e, false end,
  Str = math_str,
}
