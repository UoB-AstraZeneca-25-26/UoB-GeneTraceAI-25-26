import React, { useRef, useState, useEffect } from 'react';
import { AssistantContext, ChatMessage } from '../../types';
import { sendMessage as defaultSend, SendMessage, SUGGESTED_PROMPTS } from '../../lib/assistant';
import { Sparkles, Send, Loader2, User } from 'lucide-react';

interface AssistantPanelProps {
  context: AssistantContext | null;
  send?: SendMessage; // override to plug in your own agent
}

export const AssistantPanel: React.FC<AssistantPanelProps> = ({ context, send = defaultSend }) => {
  const [messages, setMessages] = useState<ChatMessage[]>([]);
  const [input, setInput] = useState('');
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const scrollRef = useRef<HTMLDivElement>(null);

  useEffect(() => {
    scrollRef.current?.scrollTo({ top: scrollRef.current.scrollHeight, behavior: 'smooth' });
  }, [messages, loading]);

  const submit = async (text: string) => {
    const content = text.trim();
    if (!content || loading) return;
    const next: ChatMessage[] = [...messages, { role: 'user', content }];
    setMessages(next);
    setInput('');
    setLoading(true);
    setError(null);
    try {
      const reply = await send(next, context);
      setMessages([...next, { role: 'assistant', content: reply }]);
    } catch (e) {
      setError(e instanceof Error ? e.message : 'The assistant failed to respond.');
    } finally {
      setLoading(false);
    }
  };

  const hasContext = context && context.topLines.length > 0;

  return (
    <div className="max-w-3xl mx-auto h-full flex flex-col">
      <div className="mb-4">
        <h2 className="text-2xl font-bold text-slate-900 font-display tracking-tight flex items-center gap-2">
          <Sparkles className="w-5 h-5 text-indigo-600" /> Assistant
        </h2>
        <p className="text-sm text-slate-500 mt-1">
          Ask how a ranking was produced, or anything else about your query.
        </p>
      </div>

      {/* Context chip */}
      {hasContext && (
        <div className="mb-3 inline-flex items-center gap-2 self-start bg-indigo-50 border border-indigo-200 text-indigo-800 rounded-full px-3 py-1 text-xs">
          <span className="font-semibold">Context:</span>
          <span className="font-mono">{context!.genes.join(', ')}</span>
          {context!.mode && <span className="text-indigo-400">· {context!.mode}</span>}
        </div>
      )}

      {/* Messages */}
      <div
        ref={scrollRef}
        className="flex-1 overflow-y-auto bg-white rounded-xl border border-slate-200 shadow-sm p-4 space-y-4 min-h-[320px]"
      >
        {messages.length === 0 && (
          <div className="h-full flex flex-col items-center justify-center text-center text-slate-400 py-10">
            <Sparkles className="w-8 h-8 mb-3 text-slate-300" />
            <p className="text-sm max-w-sm">
              {hasContext
                ? `I can see your current ${context!.mode ?? ''} query for ${context!.genes.join(', ')}. Ask me to explain the ranking, or anything else.`
                : 'Run a query first and I can explain the ranking — or ask me anything now.'}
            </p>
          </div>
        )}

        {messages.map((m, i) => (
          <div key={i} className={`flex gap-3 ${m.role === 'user' ? 'flex-row-reverse' : ''}`}>
            <div
              className={`shrink-0 w-8 h-8 rounded-lg flex items-center justify-center ${
                m.role === 'user' ? 'bg-slate-200 text-slate-600' : 'bg-indigo-600 text-white'
              }`}
            >
              {m.role === 'user' ? <User className="w-4 h-4" /> : <Sparkles className="w-4 h-4" />}
            </div>
            <div
              className={`rounded-xl px-4 py-2.5 text-sm max-w-[80%] whitespace-pre-wrap leading-relaxed ${
                m.role === 'user'
                  ? 'bg-slate-100 text-slate-800'
                  : 'bg-indigo-50 text-slate-800 border border-indigo-100'
              }`}
            >
              {m.content}
            </div>
          </div>
        ))}

        {loading && (
          <div className="flex gap-3">
            <div className="shrink-0 w-8 h-8 rounded-lg bg-indigo-600 text-white flex items-center justify-center">
              <Sparkles className="w-4 h-4" />
            </div>
            <div className="rounded-xl px-4 py-2.5 bg-indigo-50 border border-indigo-100">
              <Loader2 className="w-4 h-4 animate-spin text-indigo-500" />
            </div>
          </div>
        )}

        {error && <div className="text-xs text-rose-600 bg-rose-50 border border-rose-200 rounded-lg p-2">{error}</div>}
      </div>

      {/* Suggested prompts */}
      {messages.length === 0 && (
        <div className="flex flex-wrap gap-2 mt-3">
          {SUGGESTED_PROMPTS.map((p) => (
            <button
              key={p}
              onClick={() => submit(p)}
              className="text-xs bg-white border border-slate-200 hover:border-indigo-300 hover:text-indigo-700 text-slate-600 rounded-full px-3 py-1.5 transition"
            >
              {p}
            </button>
          ))}
        </div>
      )}

      {/* Input */}
      <div className="mt-3 flex items-end gap-2">
        <textarea
          value={input}
          onChange={(e) => setInput(e.target.value)}
          onKeyDown={(e) => {
            if (e.key === 'Enter' && !e.shiftKey) {
              e.preventDefault();
              submit(input);
            }
          }}
          rows={1}
          placeholder="Ask about the ranking…"
          className="flex-1 resize-none bg-white border border-slate-300 rounded-xl px-4 py-3 text-sm focus:ring-2 focus:ring-indigo-500 focus:outline-none max-h-32"
        />
        <button
          onClick={() => submit(input)}
          disabled={loading || !input.trim()}
          className="shrink-0 bg-indigo-600 hover:bg-indigo-700 disabled:opacity-40 disabled:cursor-not-allowed text-white rounded-xl p-3 transition"
          aria-label="Send"
        >
          <Send className="w-4 h-4" />
        </button>
      </div>
    </div>
  );
};