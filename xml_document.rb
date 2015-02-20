class XmlDocument
  def initialize (indent = false ,counter=-1)
    @indent = indent
    @counter = counter
  end

  def method_missing(m, args= "", &block)
    send(m, args, &block)
  end
  def newline_adder
    if @indent == true
      "\n"
    end
  end
  def indent_adder
    if @indent == true
      "  " * @counter
    else
      ''
    end
  end
  def send(m, args= "", &block)
      @counter+=1
    if block_given?
        front = "#{indent_adder}<#{m}>#{indent_check}"
        back = "#{indent_adder}</#{m}>#{newline_adder}"
      return front + block.call + back
    else
      if args.empty?
        "#{indent_adder}<#{m}/>#{newline_adder}"
      else
        argument_key = args.keys[0]
        "#{indent_adder}<#{m} #{argument_key}='#{args[argument_key]}'/>#{newline_adder}"
      end
    end
  end
end
